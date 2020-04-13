const STATUS_INFECTED = 'I';
const STATUS_RECOVERED = 'R';
const STATUS_HEALTHY = 'H';
const MAX_VELOCITY = 100;

// Credits to https://github.com/DbxDev/NBodySimulation for a large part of the code.

function StaticValues(canvas) {
    this.CONTEXT = canvas.getContext('2d');
    this.BORDER = 0;
    // Simulation values [0; 1]
    this.MIN_X = 0;
    this.MIN_Y = 0;
    this.MAX_X = 1;
    this.MAX_Y = 1;
    // Drawing values
    this.MIN_X_COORD = this.BORDER;
    this.MIN_Y_COORD = this.BORDER;
    this.MAX_X_COORD = canvas.width - this.BORDER;
    this.MAX_Y_COORD = canvas.height - this.BORDER;
    // Normalized ratio
    this.NORM_X_RATIO = (this.MAX_X_COORD - this.MIN_X_COORD) / (this.MAX_X - this.MIN_X)
    this.NORM_Y_RATIO = (this.MAX_Y_COORD - this.MIN_Y_COORD) / (this.MAX_Y - this.MIN_Y);

}

// Result between 0 and 1

function normalizedXDistance(distance) {
    return Math.round(distance * STATIC_VALUES.NORM_X_RATIO);
}

function normalizedYDistance(distance) {
    return Math.round(distance * STATIC_VALUES.NORM_Y_RATIO);
}

// Generate an ID
// Increments the IN every time it is called.

function UniqueIDGenerator(init) {

    var lastID=-1;

    if (init) {
        lastID = init - 1;
    }

    this.newId = function() {
        lastID++;
        return lastID;
    };
}

var IDGen = new UniqueIDGenerator();

function initHiddenCanvas(color, radius) {
    var hiddenCanvas;
    var hiddenContext;
    hiddenCanvas = document.createElement('canvas');
    hiddenCanvas.width = normalizedXDistance(2 * radius);
    hiddenCanvas.height = normalizedYDistance(2 * radius);
    hiddenContext = hiddenCanvas.getContext('2d');

    hiddenContext.fillStyle = color;
    hiddenContext.beginPath();
    hiddenContext.arc(
        normalizedXDistance(radius),
        normalizedYDistance(radius),
        normalizedXDistance(radius),
        0,
        Math.PI * 2); // x, y, radius, starting angle, ending angle [, option clockwise]
    hiddenContext.fill();
    hiddenContext.closePath();
    return hiddenCanvas;
}

// Sphere is the main class that controls the visualization, one sphere per particle.
// A sphere has a radius, mass, position (x, y), velocity (vx, vy), and forces acting on it (fx, fy).
// Also it remembers the forces of the previous step (fx_old, fx_new).

function Sphere(radius, mass, x, y, vx, vy, id) {

    // Except 0 we require an ID with a valid value
    if (id !== 0 && !id) {
        id = IDGen.newId();
    }
    this.id = id
    this.radius = radius;
    this.mass = mass;
    this.x = x;
    this.y = y;
    this.vx = vx;
    this.vy = vy;
    this.fx = 0.0;
    this.fy = 0.0;
    this.fx_old =  0.0;
    this.fy_old =  0.0;
    this.color = rgb(256, 256, 256); // Initial color, will be updated as soon as status changes
    this.hiddenCanvas = initHiddenCanvas(this.color,this.radius);
    this.status = undefined; //status initialized to undefined
}

// Utility function to change the status of a sphere.
// Possible states are hard-coded at the top of this file.
// TODO: figure out a good way to control also state transitions.
// I.e., no-one should go from dead to healthy or similar...

Sphere.prototype.setStatus = function(newstatus) {

    if (newstatus == this.status) {
        console.log("Already have this status");
        return
    }
    hiddenContext = this.hiddenCanvas.getContext('2d');
    if (newstatus == STATUS_HEALTHY) {
        this.status = STATUS_HEALTHY;
        hiddenContext.fillStyle = rgb(0, 255, 0);
    }
    else if (newstatus == STATUS_INFECTED) {
        this.status = STATUS_INFECTED;
        hiddenContext.fillStyle = rgb(255, 0, 0);
    }
    else if (newstatus == STATUS_RECOVERED) {
        this.status = STATUS_RECOVERED;
        hiddenContext.fillStyle = rgb(0, 0, 255);
    }
    else {
        throw new Error("Unknown status "+newstatus + ".");
    }
    // Update the colour in the animation:
    hiddenContext.fill();

};

Sphere.prototype.toString = function() {
        return "{Sphere#" + this.id
            + " R=" + this.radius
            + " m=" + this.mass
            + " (x,y)=" + this.x + "," + this.y
            + " (vx,vy)=" + this.vx + "," + this.vy+"}";
};

Sphere.prototype.Draw = function(context) {
    context.drawImage(this.hiddenCanvas,
            normalizedXDistance(this.x) - normalizedXDistance(this.radius),
            normalizedYDistance(this.y) - normalizedYDistance(this.radius));
};

// Propagating based on so-called Velocity Verlet integration.
// https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
// x(t + dt) = x(t) + v(t) * dt + 0.5 * F(t) / m * dt^2

Sphere.prototype.propagatePosVerlet = function(dt) {
    this.x = (this.x + this.vx * dt + 0.5 * this.fx / this.mass * dt**2) % 1.0;
    this.y = (this.y + this.vy * dt + 0.5 * this.fy / this.mass * dt**2) % 1.0;
}

// Propagating velocities based on Velocity Verlet integration.
// https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
// v(t + dt) = v(t) + dt * [a(t) + a(t + dt)] / 2

Sphere.prototype.propagateVelVerlet = function(dt) {
    this.vx = this.vx + dt * (this.fx_old + this.fx) / (2 * this.mass);
    this.vy = this.vy + dt * (this.fy_old + this.fy) / (2 * this.mass);

    // Checks on the velocity to have a bound on them.
    // This is very dirty and breaks energy conservation, but needed to make sure the
    // cell does not explode when the time step is too large and sliders are moved very fast.
    // TODO: find better way to deal with this.
    if (this.vx > MAX_VELOCITY) this.vx = MAX_VELOCITY;
    if (this.vx < -MAX_VELOCITY) this.vx = -MAX_VELOCITY;
    if (this.vy > MAX_VELOCITY) this.vy = MAX_VELOCITY;
    if (this.vy < -MAX_VELOCITY) this.vy = -MAX_VELOCITY;

}

function rgb(r,g,b) {
    return "rgb(" + r + "," + g + "," + b + ")";
}
