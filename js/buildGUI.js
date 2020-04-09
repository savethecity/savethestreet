var freqHandwashingSlider = document.getElementById("frequencyHandWashing");
var socialDistancingSlider = document.getElementById("socialDistancing");
var economicActivitySlider = document.getElementById("economicActivity");


function sleep(ms) {
  // need this function to pause the simulation for ms times at every step.
  // TODO: is this the best way to do it?
  return new Promise(resolve => setTimeout(resolve, ms));
}

function shuffleArray( array ) {
    //shuffling the elements of the array.
    var size = array.length;
    var rand, tmp;
    for (var i=0 ; i<size ; i++ ) {
        rand = Math.floor(Math.random()*(size-i))
        tmp = array[i];
        array[i] = array[rand];
        array[rand] = tmp;
    }
    return array;
}

function generateNSpheres(N , R) {
    /* Generates N spheres of radius R.
    When initializing, it is not smart to initialize positions and velocities randomly!
    Particles might be too close to each other, and this will lead to very large accelerations.
    current Algorithm: Generate a uniform grid, that is based on the radius, so that spheres of radius R
    can be put on neighboring gridpoints with at least 2R distance, plus some margin.
    I populate the gridpoints randomly with N spheres. This leads to a well ordered initial state,
    which looks strange. TODO: find better way to populate in the beginning */

    var r_margin= R * 1.2 
    var D = 2 * R ;
    var d_margin = 2*R * 1.5 ; //50% margin between 2 objects
    var elem_on_one_line = parseInt(1/d_margin);
    var maxIndex = elem_on_one_line **2;
    const spheres = new Array();

    if (maxIndex < N) throw new Error("Impossible situation, too many or too big spheres.");

    // Let's keep a bound on the computational cost of this:
    if (maxIndex>1000000) throw new Error("Radius is too small");

    occupied = [];
    for (var i=0 ; i<N ; i++) occupied[i]=true; // reservation of N values for N spheres
    // fill the rest of the array with false.
    for (i= N ; i< maxIndex ; i++) occupied[i] = false;
    // Shuffle array so that occupied[i]=true is distributed randomly
    shuffleArray(occupied);
    var count = 0;
    for (var i=0 ; i<maxIndex ; i++) {
        if (occupied[i]) {
            // Random velocity. TODO: take velocity from Boltzmann distribution at set temperature!
            vx= (1-2*Math.random())*0.3;
            vy=(1-2*Math.random())*0.3;
            x = r_margin + d_margin*Math.floor(i/elem_on_one_line);
            y = r_margin+d_margin*(i%elem_on_one_line);
            sphere = new Sphere(R, 1, x, y, vx, vy, count);
            sphere.setStatus(STATUS_HEALTHY);
            spheres[count] = sphere;
            count++;
        }
    }
    // Setting the first sphere as infected
    spheres[N/2].setStatus(STATUS_INFECTED)
    return spheres
}

function clearCanvas(){
    STATIC_VALUES.CONTEXT.clearRect(STATIC_VALUES.MIN_X_COORD,
        STATIC_VALUES.MIN_Y_COORD, STATIC_VALUES.MAX_X_COORD, STATIC_VALUES.MAX_Y_COORD);
}

function drawSpheres(spheres){
    for (var i = 0; i < spheres.length; i++) {
        spheres[i].Draw(STATIC_VALUES.CONTEXT)
    }
}

function calculateForces(spheres, epsilon, rm, infection_threshold){
    /*
    This function loops over all pairs of particles and calculates the forces on particles based
    on the LJ parameters epsilon and rm. TODO: more efficient than looping over all pairs? Neighbor lists
    The infection is propagated in this function since the distance is calculated, and the infection rate (for now)
    is distance dependant. TODO: improve, it's not great that this is shuffled into this function here
    */
    // Variables for calculation of distance
    var dist_x, dist_y, dist, rm_over_dist;
    // variables for enumeration
    var i, j;
    // force value
    var force;
    // Setting energies and forces to 0 before start of summation over all pairs:
    // ener refers to the total potential energy
    var ener=0.0;
    for (i=0; i<spheres.length; i++){
        //  copying forces into another array because they are needed for velocity propagation rusing
        // velocity verlet.
        spheres[i].old_fx = spheres[i].fx;
        spheres[i].old_fy = spheres[i].fy;
        spheres[i].fx = 0.0;
        spheres[i].fy = 0.0;
    }
    for (i=0; i<spheres.length; i++){
        for (j=i+1; j<spheres.length; j++){
            dist_x = spheres[i].x - spheres[j].x;
            dist_y = spheres[i].y - spheres[j].y;
            dist = Math.sqrt(dist_x**2 + dist_y**2);
            rm_over_dist = (rm / dist)
            ener += epsilon*( (rm / dist)**12 - 2* (rm / dist)**6);
            force = - 12*epsilon * ( rm**12 / dist**14  - rm**6  / dist**8)
            //~ force = -12 * A / dist**14 * dist_x + 6 *B / dist**8 * dist_x;
            spheres[i].fx -= force*dist_x
            spheres[j].fx += force*dist_x
            //~ force = -12 * A / dist**14 * dist_y + 6 *B / dist**8 * dist_y;
            spheres[i].fy -= force*dist_y
            spheres[j].fy += force*dist_y
            
            /* Misusing the fact that distance has been calculated in order to process infections,
            since that is distance dependant in this model. p_infection ~ e^{-dist/infection_threshold}
            */
            if ( dist < 3*infection_threshold) {
                if ((spheres[i].status == STATUS_INFECTED) && (spheres[j].status == STATUS_HEALTHY)) {
                    if (Math.random() < Math.exp(- dist / infection_threshold)){
                            spheres[j].setStatus(STATUS_INFECTED)
                            console.log('Infecting '+j)
                    }
                } else if ((spheres[j].status == STATUS_INFECTED) && (spheres[i].status == STATUS_HEALTHY)) {
                    if (Math.random() < Math.exp(- dist / infection_threshold)){
                            spheres[i].setStatus(STATUS_INFECTED)
                            console.log('Infecting '+i)
                    }
                }
            }
                    
        }
    }
    return ener
}

function propagatePositionsVerlet(spheres, dt){
    var i;
    var sphere;
    for (i = 0; i < spheres.length; i++) {
        sphere = spheres[i]
        sphere.propagatePosVerlet(dt)
    }
}

function propagateVelocitiesVerlet(spheres, dt){
    var i;
    var sphere;
    for (i = 0; i < spheres.length; i++) {
        sphere = spheres[i]
        sphere.propagateVelVerlet(dt)
        // Bouncing off the walls is done here:
        if ( (sphere.x < sphere.radius) && (sphere.vx < 0)) {
            sphere.vx = -sphere.vx
        }
        else if ( (sphere.x > (1-sphere.radius)) && (sphere.vx > 0)) {
            sphere.vx = -sphere.vx
        }
        if ( (sphere.y < sphere.radius) && (sphere.vy < 0)) {
            sphere.vy = -sphere.vy
        }
        else if ((sphere.y > (1-sphere.radius)) && (sphere.vy> 0)) {
            sphere.vy = -sphere.vy
        }
    }
}


function getKineticEnergy(spheres){
    var ekin = 0.0;
    for (var i=0; i<spheres.length; i++){
        ekin += 0.5 * spheres[i].mass * (spheres[i].vx**2 +spheres[i].vy**2)
    }
    return ekin
}

function applyThermostat(spheres, ekin, target){
    /* Very basic velocity rescaling. Rescales all velocities uniformly so that new kinetic energy
    is at a target value.
    
     */
    var aux = Math.sqrt(target/ekin)
    for (var i = 0; i < spheres.length; i++) {
        spheres[i].vx *= aux
        spheres[i].vy *= aux
    }
}

function timelineInfected(spheres, p_recovery, p_fatality){
    for (var i = 0; i < spheres.length; i++) {
        if (spheres[i].status == STATUS_INFECTED){
            if (Math.random() < p_recovery){
                spheres[i].setStatus(STATUS_RECOVERED)
            }
        }
    }
}

var STATIC_VALUES = null;
window.onload = async function()
{
    var canvas = document.getElementById('my_canvas');
        if(!canvas)
        {
            alert("Canvas not found.");
            return;
        }
    var context = canvas.getContext('2d');
        if(!context)
        {
            alert("Canvas context not found.");
            return;
        }
    var epsilon = 1e-13;
    var rm = 0.2;
    var infection_threshold, temperature;
    var count_sick;
    const dt = 0.0005;
    var i, j, epot, etot;

    STATIC_VALUES = new StaticValues(canvas);

    spheres = generateNSpheres(600, 0.01);
    // drawing the spheres immeiately before doing anything else
    drawSpheres(spheres)

    // Calculating forces and potential energy in this first step:
    epot = calculateForces(spheres, epsilon, rm)
    //~ ekin = getKineticEnergy(spheres) 
    //~ etot = (ekin+epot)         
    istep = 0;
    while (true) {
        // Getting all the slide values here:
        // This is where the "physics" (or rather psychology) needs to come in!
        // How does social distancing affect infection rates?
        rm = 0.1 + parseInt(socialDistancingSlider.value)/400; // social distancing slider
        // How does the handwashing frequency affect infection rates?
        infection_threshold = 1.0/parseInt(freqHandwashingSlider.value)
        // Economic activity gives the "temperature" of the simulation!.
        temperature = parseInt(economicActivitySlider.value)

        // propagating spheres based on velocity verlet scheme:
        propagatePositionsVerlet(spheres, dt)
        // Now I need the forces on the new positions x(t+dt):
        epot = calculateForces(spheres, epsilon, rm, infection_threshold)
        // propagating velocities based on old velocities&forces and new forces:
        propagateVelocitiesVerlet(spheres, dt)
        // get the kinetic energy from v(t+dt)
        ekin = getKineticEnergy(spheres) 
        
        // applying the thermostat which is parametrized by the economic activity slider:
        //~ if ( (istep % 10) == 0 ){
        // Thermostats needs to be apply every now and then
        applyThermostat(spheres, ekin, temperature)

        // etot(t+dt) = ekin(t+dt) + epot(t+dt)
        etot = (ekin+epot)
        console.log(istep+': '+ etot + ' ' + ekin + '   '+epot+'   ' +rm)
        timelineInfected(spheres, 0.001, 0.01)
        await sleep(20)
        clearCanvas()
        drawSpheres(spheres)
        istep++;
        if ((istep % 100) == 0) {
            console.log('Checking number of sick')
            count_sick = 0;
            for (var i=0; i<spheres.length; i++){
                if (spheres[i].status ==  STATUS_INFECTED) count_sick++;
            }
            console.log(count_sick+' sick individuals')
            if (count_sick == 0) break // stopping if there's noone infected.
        }
    }
}
