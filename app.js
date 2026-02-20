// Twin Timer with robust gravity removal + rotation gating
// Key changes:
// - Use Madgwick IMU (gyro+acc) or AbsoluteOrientationSensor for quaternion.
// - Transform raw acceleration (including gravity) to world frame,
//   low-pass filter to estimate gravity (~ [0,0,-g]), subtract to get linear a_world.
// - Gate integration during rotation: if |omega| > thresh and |a| small -> skip integrate, decay v.
// - Adaptive drift: tau_eff = tau / (1 + k*(gamma - 1)).

const els = {
  timer: document.getElementById('timer'),
  startStop: document.getElementById('startStop'),
  reset: document.getElementById('reset'),
  openSettings: document.getElementById('openSettings'),
  closeSettings: document.getElementById('closeSettings'),
  closeSettingsFooter: document.getElementById('closeSettingsFooter'),
  modalBackdrop: document.getElementById('modalBackdrop'),
  settings: document.getElementById('settings'),
  accel: document.getElementById('accel'),
  vraw: document.getElementById('vraw'),
  veff: document.getElementById('veff'),
  gamma: document.getElementById('gamma'),
  oriStatus: document.getElementById('oriStatus'),
  cval: document.getElementById('cval'),
  tau: document.getElementById('tau'),
  kboost: document.getElementById('kboost'),
};

// ---- Config getters ----
function getC() {
  const c = parseFloat(els.cval.value);
  return Number.isFinite(c) && c > 0 ? c : 1.0;
}
function getTau() {
  const t = parseFloat(els.tau.value);
  return Number.isFinite(t) && t > 0 ? t : 3.0;
}
function getK() {
  const k = parseFloat(els.kboost.value);
  return Number.isFinite(k) && k >= 0 ? k : 0.5;
}

// ---- Math & helpers ----
const G = 9.80665;                 // m/s^2
const ACC_GATE = 0.05;             // m/s^2, noise gate
const OMEGA_THRESH = 0.35;         // rad/s, rotation gating threshold
const GRAV_TAU = 0.7;              // s, gravity LPF time-constant in world frame

const norm3 = (x, y, z) => Math.hypot(x, y, z);
const addScaled = (v, a, dt) => ({ x: v.x + a.x * dt, y: v.y + a.y * dt, z: v.z + a.z * dt });

function formatTau(sec) {
  const ms = Math.floor((sec % 1) * 1000);
  const total = Math.floor(sec);
  const s = total % 60;
  const m = Math.floor((total / 60) % 60);
  const h = Math.floor(total / 3600);
  const pad = (n, w = 2) => n.toString().padStart(w, '0');
  return `${pad(h)}:${pad(m)}:${pad(s)}.${pad(ms, 3)}`;
}
function limitSpeed(vmag, c) { return c * Math.tanh(vmag / c); }
function gammaFromV(vmag, c) {
  const beta2 = (vmag * vmag) / (c * c);
  if (beta2 <= 0) return 1;
  if (beta2 >= 1) return 1e12;
  return 1 / Math.sqrt(1 - beta2);
}

// ---- Quaternion ops: q = [w,x,y,z] ----
function qNormalize(q) {
  const n = Math.hypot(q[0], q[1], q[2], q[3]) || 1;
  return [q[0]/n, q[1]/n, q[2]/n, q[3]/n];
}
function qConj(q) { return [q[0], -q[1], -q[2], -q[3]]; }
function qMul(a, b) {
  return [
    a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3],
    a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2],
    a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1],
    a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0],
  ];
}
function qRotateVec(q, v) {
  const vq = [0, v.x, v.y, v.z];
  const t = qMul(q, vq);
  const r = qMul(t, qConj(q));
  return { x: r[1], y: r[2], z: r[3] };
}

// ---- Madgwick IMU (gyro+acc only) for orientation (world->device) ----
const MadgwickIMU = (function() {
  // beta trades gyro drift vs accel trust. 0.02..0.1 typical. Lower = smoother.
  const beta = 0.05;
  let qwd = [1,0,0,0]; // world->device
  function reset() { qwd = [1,0,0,0]; }
  function getQwd() { return qwd; }
  function update(gx, gy, gz, ax, ay, az, dt) {
    // Normalize accel
    const aNorm = Math.hypot(ax, ay, az);
    if (aNorm === 0) return;
    ax /= aNorm; ay /= aNorm; az /= aNorm;

    let q1=qwd[0], q2=qwd[1], q3=qwd[2], q4=qwd[3];

    // Objective: gravity should align with device-frame -z axis
    // Compute gradient descent step (see Madgwick 2010, IMU update)
    const _2q1 = 2*q1, _2q2 = 2*q2, _2q3 = 2*q3, _2q4 = 2*q4;
    const _4q1 = 4*q1, _4q2 = 4*q2, _4q3 = 4*q3;
    const q1q1 = q1*q1, q2q2 = q2*q2, q3q3 = q3*q3, q4q4 = q4*q4;

    const f1 = _2q2*q4 - _2q1*q3 - ax;
    const f2 = _2q1*q2 + _2q3*q4 - ay;
    const f3 = 1 - 2*(q2q2 + q3q3) - az;

    const s1 = -_2q3*f1 + _2q2*f2;
    const s2 =  _2q4*f1 + _2q1*f2 - _4q2*f3;
    const s3 = -_4q3*f3 - _2q1*f1 + _2q4*f2;
    const s4 =  _2q2*f1 + _2q3*f2;

    const recip = 1 / Math.hypot(s1, s2, s3, s4 || 1);
    const s1n = s1 * recip, s2n = s2 * recip, s3n = s3 * recip, s4n = s4 * recip;

    // Quaternion rate from gyro + feedback
    const qDot1 = 0.5 * (-q2*gx - q3*gy - q4*gz) - beta * s1n;
    const qDot2 = 0.5 * ( q1*gx + q3*gz - q4*gy) - beta * s2n;
    const qDot3 = 0.5 * ( q1*gy - q2*gz + q4*gx) - beta * s3n;
    const qDot4 = 0.5 * ( q1*gz + q2*gy - q3*gx) - beta * s4n;

    // Integrate and normalize
    q1 += qDot1 * dt;
    q2 += qDot2 * dt;
    q3 += qDot3 * dt;
    q4 += qDot4 * dt;
    qwd = qNormalize([q1, q2, q3, q4]);
  }
  return { update, reset, getQwd };
})();

// ---- Sensor state ----
let running = false;
let lastFrameT = null;

let linAccSensor = null;   // LinearAccelerationSensor
let accSensor = null;      // Accelerometer (includes gravity)
let gyroSensor = null;     // Gyroscope
let absOriSensor = null;   // AbsoluteOrientationSensor

let hasAbsQ = false;
let qdw_abs = [1,0,0,0];   // device->world from AbsoluteOrientationSensor
let lastGyroTS = null;
let omega = { x:0, y:0, z:0 }; // rad/s

// Gravity estimator in world frame (low-pass)
let gWorld = { x: 0, y: 0, z: -G }; // start with -z

// World-frame velocity
let vWorld = { x: 0, y: 0, z: 0 };

// Proper time
let tauVal = 0;

// ---- Start sensors (Generic Sensor API preferred) ----
function startGenericSensors() {
  try {
    if ('AbsoluteOrientationSensor' in window) {
      absOriSensor = new AbsoluteOrientationSensor({ frequency: 60, referenceFrame: 'device' });
      absOriSensor.addEventListener('reading', () => {
        const q = absOriSensor.quaternion;
        if (q) {
          qdw_abs = [q[0], q[1], q[2], q[3]];
          hasAbsQ = true;
          els.oriStatus.textContent = 'absolute quaternion';
        }
      });
      absOriSensor.start();
    }

    if ('LinearAccelerationSensor' in window) {
      linAccSensor = new LinearAccelerationSensor({ frequency: 60 });
      linAccSensor.addEventListener('reading', () => { /* read in loop */ });
      linAccSensor.start();
    } else if ('Accelerometer' in window) {
      accSensor = new Accelerometer({ frequency: 60 });
      accSensor.addEventListener('reading', () => { /* read in loop */ });
      accSensor.start();
    }

    if ('Gyroscope' in window) {
      gyroSensor = new Gyroscope({ frequency: 60 });
      gyroSensor.addEventListener('reading', () => {
        const t = gyroSensor.timestamp;
        omega = { x: gyroSensor.x || 0, y: gyroSensor.y || 0, z: gyroSensor.z || 0 };
        if (!hasAbsQ) {
          if (accSensor || linAccSensor) {
            if (lastGyroTS != null) {
              const dt = Math.max(0.0001, (t - lastGyroTS) / 1000);
              // For IMU fusion, use accel including gravity (best to use accSensor) OR linear if present
              const aDev = getAccelIncludingGravityDev();
              MadgwickIMU.update(omega.x, omega.y, omega.z, aDev.x, aDev.y, aDev.z, dt);
              els.oriStatus.textContent = 'IMU fusion';
            }
          }
        }
        lastGyroTS = t;
      });
      gyroSensor.start();
    }
    return true;
  } catch {
    return false;
  }
}

// ---- Legacy fallbacks (DeviceMotion/DeviceOrientation) ----
let dmAccelInc = { x:0,y:0,z:0 };     // including gravity
let dmLinAcc = { x:0,y:0,z:0 };       // without gravity if provided
let dmOmega = { x:0,y:0,z:0 };
function startLegacySensors() {
  let ok = false;
  if ('DeviceMotionEvent' in window) {
    window.addEventListener('devicemotion', (e) => {
      const ag = e.accelerationIncludingGravity;
      if (ag && ag.x != null) dmAccelInc = { x: ag.x||0, y: ag.y||0, z: ag.z||0 };
      const a = e.acceleration;
      if (a && a.x != null) dmLinAcc = { x: a.x||0, y: a.y||0, z: a.z||0 };
      if (e.rotationRate) {
        // DeviceMotion rotationRate: alpha (z), beta (x), gamma (y) in deg/s
        dmOmega = {
          x: (e.rotationRate.beta || 0) * Math.PI/180,
          y: (e.rotationRate.gamma || 0) * Math.PI/180,
          z: (e.rotationRate.alpha || 0) * Math.PI/180
        };
      }
      // IMU update driven by display loop timing in our case
    }, { passive: true });
    ok = true;
  }
  if ('DeviceOrientationEvent' in window) {
    // We avoid Euler here; rely on IMU to reduce gravity leakage from yaw noise.
    ok = true;
  }
  return ok;
}

// ---- Sensor reading helpers ----
function getAccelIncludingGravityDev() {
  if (accSensor && accSensor.x != null) {
    return { x: accSensor.x || 0, y: accSensor.y || 0, z: accSensor.z || 0 };
  }
  // Legacy
  return dmAccelInc;
}
function getLinearAccelDev() {
  if (linAccSensor && linAccSensor.x != null) {
    return { x: linAccSensor.x || 0, y: linAccSensor.y || 0, z: linAccSensor.z || 0 };
  }
  // Legacy linear accel if present:
  if (dmLinAcc) return dmLinAcc;
  return null;
}
function getOmega() {
  if (gyroSensor && gyroSensor.x != null) return omega;
  return dmOmega;
}
function getQdw() {
  if (hasAbsQ) return qdw_abs;
  // IMU gives world->device; invert to get device->world
  return qConj(MadgwickIMU.getQwd());
}

// ---- Main loop ----
let lastLoopT = null;
function loop(t) {
  if (!running) return;
  if (lastLoopT == null) lastLoopT = t;
  const dt = Math.max(0.0005, (t - lastLoopT) / 1000);
  lastLoopT = t;

  const qdw = getQdw();
  const omg = getOmega();
  const omegaMag = norm3(omg.x, omg.y, omg.z);

  // Build world-frame linear acceleration
  let aWorld = { x:0, y:0, z:0 };

  const aLinDev = getLinearAccelDev();
  if (aLinDev) {
    // If we have linear accel directly, just rotate to world
    aWorld = qRotateVec(qdw, aLinDev);
  } else {
    // Otherwise, use accel including gravity -> world, subtract estimated gravity
    const aIncDev = getAccelIncludingGravityDev();
    const aIncWorld = qRotateVec(qdw, aIncDev);

    // Gravity estimator (world): gWorld <- alpha*gWorld + (1-alpha)*aIncWorld, then normalize to |g|
    const alpha = Math.exp(-dt / GRAV_TAU);
    gWorld.x = alpha * gWorld.x + (1 - alpha) * aIncWorld.x;
    gWorld.y = alpha * gWorld.y + (1 - alpha) * aIncWorld.y;
    gWorld.z = alpha * gWorld.z + (1 - alpha) * aIncWorld.z;
    // Enforce magnitude ~ G to suppress drift
    const gmag = norm3(gWorld.x, gWorld.y, gWorld.z) || G;
    const scale = G / gmag;
    gWorld.x *= scale; gWorld.y *= scale; gWorld.z *= scale;

    aWorld = { x: aIncWorld.x - gWorld.x, y: aIncWorld.y - gWorld.y, z: aIncWorld.z - gWorld.z };
  }

  // Noise gate in world frame
  if (Math.abs(aWorld.x) < ACC_GATE) aWorld.x = 0;
  if (Math.abs(aWorld.y) < ACC_GATE) aWorld.y = 0;
  if (Math.abs(aWorld.z) < ACC_GATE) aWorld.z = 0;

  // Rotation gating: if rotating but not accelerating, don't integrate; decay velocity instead
  const aMag = norm3(aWorld.x, aWorld.y, aWorld.z);
  const rotatingButStill = (omegaMag > OMEGA_THRESH) && (aMag < 0.2);
  if (!rotatingButStill) {
    // Integrate to velocity
    vWorld = addScaled(vWorld, aWorld, dt);
  }

  // Relativistic speed handling
  const vmag = norm3(vWorld.x, vWorld.y, vWorld.z);
  const c = getC();
  const vEff = limitSpeed(vmag, c);
  const gamma = gammaFromV(vEff, c);

  // Drift: stronger when still and/or rotating-only; adaptive with gamma
  const tau = getTau();
  const tauEff = tau / (1 + getK() * (gamma - 1));
  let decay = Math.exp(-dt / tauEff);
  if (rotatingButStill) {
    // Extra decay during rotation-only to shed spurious v quickly
    decay = Math.exp(-dt / Math.max(0.3 * tauEff, 0.05));
  }
  // Also apply decay when acceleration is tiny (stillness)
  if (aMag === 0 || aMag < ACC_GATE * 1.5) {
    vWorld.x *= decay; vWorld.y *= decay; vWorld.z *= decay;
  }

  // Proper time update: Δτ = Δt / γ
  tauVal += dt / gamma;

  // UI
  els.timer.textContent = formatTau(tauVal);
  els.accel.textContent = `${aWorld.x.toFixed(3)}, ${aWorld.y.toFixed(3)}, ${aWorld.z.toFixed(3)}`;
  els.vraw.textContent = vmag.toFixed(3);
  els.veff.textContent = vEff.toFixed(3);
  els.gamma.textContent = gamma.toFixed(6);

  requestAnimationFrame(loop);
}

// ---- Controls & modal ----
function openSettings() {
  els.modalBackdrop.hidden = false;
  if (typeof els.settings.showModal === 'function') els.settings.showModal();
  else els.settings.style.display = 'block';
}
function closeSettings() {
  els.modalBackdrop.hidden = true;
  if (typeof els.settings.close === 'function') els.settings.close();
  else els.settings.style.display = 'none';
}

els.openSettings.addEventListener('click', openSettings);
els.closeSettings.addEventListener('click', closeSettings);
els.closeSettingsFooter.addEventListener('click', closeSettings);
els.modalBackdrop.addEventListener('click', closeSettings);

els.startStop.addEventListener('click', () => {
  if (!running) {
    const started = startGenericSensors() || startLegacySensors();
    els.oriStatus.textContent = started ? 'starting…' : 'no sensors';
    running = true;
    lastLoopT = null;
    requestAnimationFrame(loop);
    els.startStop.textContent = 'Pause';
  } else {
    running = false;
    els.startStop.textContent = 'Start';
  }
});

els.reset.addEventListener('click', () => {
  tauVal = 0;
  vWorld = { x: 0, y: 0, z: 0 };
  gWorld = { x: 0, y: 0, z: -G };
  els.timer.textContent = formatTau(0);
  els.vraw.textContent = '0.00';
  els.veff.textContent = '0.00';
  els.gamma.textContent = '1.000000';
});

document.addEventListener('visibilitychange', () => {
  if (document.hidden && running) {
    running = false;
    els.startStop.textContent = 'Start';
  }
});
