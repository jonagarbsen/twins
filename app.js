// Twin Timer with quaternion fusion (Madgwick) and adaptive drift
// - Fuses gyro + accel + magnetometer into a quaternion (world<->device).
// - Uses LinearAccelerationSensor if available; else subtracts gravity using quaternion.
// - Rotates linear acceleration to a fixed world frame before integrating.
// - Adaptive drift: tau_eff = tau / (1 + k*(gamma - 1)) for faster decay at high gamma.

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

// --------- Config getters ----------
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

// --------- Basic math ----------
const G = 9.80665; // m/s^2
const GATE = 0.05; // accel noise gate in device frame
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

// --------- Quaternion helpers ---------
// Quaternion q = [w, x, y, z]
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
  // v_world = q ⊗ v ⊗ q*
  const vq = [0, v.x, v.y, v.z];
  const t = qMul(q, vq);
  const r = qMul(t, qConj(q));
  return { x: r[1], y: r[2], z: r[3] };
}

// --------- Madgwick filter (world->device quaternion) ---------
// We track qwd: rotates world frame into device frame.
// To get device->world for rotating vectors, we use qdw = qConj(qwd).
const Madgwick = (function() {
  // beta controls gyro drift vs accel/mag trust; 0.05..0.2 typical
  const beta = 0.1;
  let qwd = [1,0,0,0]; // world->device
  function reset() { qwd = [1,0,0,0]; }
  function getQwd() { return qwd; }
  function update(gx, gy, gz, ax, ay, az, mx, my, mz, dt) {
    // Convert to unit vectors where appropriate
    const aNorm = Math.hypot(ax, ay, az);
    if (aNorm === 0) return;
    ax /= aNorm; ay /= aNorm; az /= aNorm;

    const mNorm = Math.hypot(mx, my, mz);
    if (mNorm !== 0) { mx /= mNorm; my /= mNorm; mz /= mNorm; }

    let q1=qwd[0], q2=qwd[1], q3=qwd[2], q4=qwd[3];

    // Reference: Madgwick 2010 report, full AHRS update
    // Auxiliary variables to avoid repeated arithmetic
    const _2q1mx = 2*q1*mx, _2q1my = 2*q1*my, _2q1mz = 2*q1*mz;
    const _2q2mx = 2*q2*mx;
    const _2q1 = 2*q1, _2q2 = 2*q2, _2q3 = 2*q3, _2q4 = 2*q4;
    const q1q1 = q1*q1, q2q2 = q2*q2, q3q3 = q3*q3, q4q4 = q4*q4;

    // Compute reference direction of Earth's magnetic field
    const hx = mx*q1q1 - _2q1my*q4 + _2q1mz*q3 + mx*q2q2 + _2q2*my*q3 + _2q2*mz*q4 - mx*q3q3 - mx*q4q4;
    const hy = _2q1mx*q4 + my*q1q1 - _2q1mz*q2 + _2q2mx*q3 - my*q2q2 + my*q3q3 + _2q3*mz*q4 - my*q4q4;
    const _2bx = Math.sqrt(hx*hx + hy*hy);
    const _2bz = -_2q1mx*q3 + _2q1my*q2 + mz*q1q1 + _2q2mx*q4 - mz*q2q2 + _2q3*my*q4 - mz*q3q3 + mz*q4q4;

    const _4bx = 2*_2bx;
    const _4bz = 2*_2bz;

    // Gradient descent algorithm corrective step
    const s1 = -_2q3*(2*(q2*q4 - q1*q3) - ax) + _2q2*(2*(q1*q2 + q3*q4) - ay)
               - _2bz*q3*(_2bx*(0.5 - q3q3 - q4q4) + _2bz*(q2*q4 - q1*q3) - mx)
               + (-_2bx*q4 + _2bz*q2)*(_2bx*(q2*q3 - q1*q4) + _2bz*(q1*q2 + q3*q4) - my)
               + _2bx*q3*(_2bx*(q1*q3 + q2*q4) + _2bz*(0.5 - q2q2 - q3q3) - mz);
    const s2 =  _2q4*(2*(q2*q4 - q1*q3) - ax) + _2q1*(2*(q1*q2 + q3*q4) - ay)
               - 4*q2*(1 - 2*(q2q2 + q3q3) - az)
               + _2bz*q4*(_2bx*(0.5 - q3q3 - q4q4) + _2bz*(q2*q4 - q1*q3) - mx)
               + (_2bx*q3 + _2bz*q1)*(_2bx*(q2*q3 - q1*q4) + _2bz*(q1*q2 + q3*q4) - my)
               + (_2bx*q4 - _4bz*q2)*(_2bx*(q1*q3 + q2*q4) + _2bz*(0.5 - q2q2 - q3q3) - mz);
    const s3 = -_2q1*(2*(q2*q4 - q1*q3) - ax) + _2q4*(2*(q1*q2 + q3*q4) - ay)
               - 4*q3*(1 - 2*(q2q2 + q3q3) - az)
               + (-_4bx*q3 - _2bz*q1)*(_2bx*(0.5 - q3q3 - q4q4) + _2bz*(q2*q4 - q1*q3) - mx)
               + (_2bx*q2 + _2bz*q4)*(_2bx*(q2*q3 - q1*q4) + _2bz*(q1*q2 + q3*q4) - my)
               + (_2bx*q1 - _4bz*q3)*(_2bx*(q1*q3 + q2*q4) + _2bz*(0.5 - q2q2 - q3q3) - mz);
    const s4 =  _2q2*(2*(q2*q4 - q1*q3) - ax) + _2q3*(2*(q1*q2 + q3*q4) - ay)
               + (-_4bx*q4 + _2bz*q2)*(_2bx*(0.5 - q3q3 - q4q4) + _2bz*(q2*q4 - q1*q3) - mx)
               + (-_2bx*q1 + _2bz*q3)*(_2bx*(q2*q3 - q1*q4) + _2bz*(q1*q2 + q3*q4) - my)
               + _2bx*q2*(_2bx*(q1*q3 + q2*q4) + _2bz*(0.5 - q2q2 - q3q3) - mz);
    const recipNorm = 1 / Math.hypot(s1, s2, s3, s4 || 1);
    const s1n = s1 * recipNorm, s2n = s2 * recipNorm, s3n = s3 * recipNorm, s4n = s4 * recipNorm;

    // Rate of change of quaternion from gyroscope
    const qDot1 = 0.5 * (-q2*gx - q3*gy - q4*gz) - beta * s1n;
    const qDot2 = 0.5 * ( q1*gx + q3*gz - q4*gy) - beta * s2n;
    const qDot3 = 0.5 * ( q1*gy - q2*gz + q4*gx) - beta * s3n;
    const qDot4 = 0.5 * ( q1*gz + q2*gy - q3*gx) - beta * s4n;

    // Integrate to yield quaternion
    q1 += qDot1 * dt;
    q2 += qDot2 * dt;
    q3 += qDot3 * dt;
    q4 += qDot4 * dt;
    qwd = qNormalize([q1, q2, q3, q4]);
  }
  return { update, reset, getQwd };
})();

// --------- Sensors and data flow ---------
let running = false;
let lastFrameT = null;

// Latest raw sensor readings
let accInc = null;    // from Accelerometer (includes gravity)
let linAcc = null;    // from LinearAccelerationSensor (without gravity)
let gyro = null;      // rad/s
let mag = null;       // in uT
let lastGyroTS = null;

// Proper time
let tauVal = 0;

// World velocity
let vWorld = { x: 0, y: 0, z: 0 };

// Orientation source flags
let usingAbsoluteOrientation = false;
let qdw_abs = [1,0,0,0]; // device->world from AbsoluteOrientationSensor

// Try Generic Sensor API first
function startGenericSensors() {
  try {
    let absOri = null, las = null, acc = null, gyr = null, magn = null;

    if ('AbsoluteOrientationSensor' in window) {
      absOri = new AbsoluteOrientationSensor({ frequency: 60, referenceFrame: 'device' });
      absOri.addEventListener('reading', () => {
        const q = absOri.quaternion; // device->world
        if (q) {
          qdw_abs = [q[0], q[1], q[2], q[3]];
          usingAbsoluteOrientation = true;
          els.oriStatus.textContent = 'absolute (quaternion)';
        }
      });
      absOri.addEventListener('error', () => { /* ignore */ });
      absOri.start();
    }

    if ('LinearAccelerationSensor' in window) {
      las = new LinearAccelerationSensor({ frequency: 60 });
      las.addEventListener('reading', () => {
        linAcc = { x: las.x || 0, y: las.y || 0, z: las.z || 0 };
      });
      las.start();
    } else if ('Accelerometer' in window) {
      acc = new Accelerometer({ frequency: 60 });
      acc.addEventListener('reading', () => {
        accInc = { x: acc.x || 0, y: acc.y || 0, z: acc.z || 0 };
      });
      acc.start();
    }

    if ('Gyroscope' in window) {
      gyr = new Gyroscope({ frequency: 60 });
      gyr.addEventListener('reading', () => {
        // Gyro rates in rad/s
        const t = gyr.timestamp;
        gyro = { x: gyr.x || 0, y: gyr.y || 0, z: gyr.z || 0, t };
        if (lastGyroTS != null) {
          const dt = Math.max(0.0001, (t - lastGyroTS) / 1000);
          // Use fusion if we have accel and mag
          const a = linAcc ?? accInc; // will normalize inside filter
          if (a && mag) {
            Madgwick.update(
              gyro.x, gyro.y, gyro.z,
              a.x, a.y, a.z,
              mag.x, mag.y, mag.z,
              dt
            );
            usingAbsoluteOrientation = false;
            els.oriStatus.textContent = 'Madgwick (fusion)';
          }
        }
        lastGyroTS = t;
      });
      gyr.start();
    }

    if ('Magnetometer' in window) {
      magn = new Magnetometer({ frequency: 30 });
      magn.addEventListener('reading', () => {
        mag = { x: magn.x || 0, y: magn.y || 0, z: magn.z || 0 };
      });
      magn.start();
    }

    return true;
  } catch (e) {
    // Some sensors may throw SecurityError if not allowed by Permissions-Policy.
    return false;
  }
}

// Fallback to DeviceMotion/DeviceOrientation if Generic Sensor not available
let dmAccel = null, dmRotRate = null, dmTS = null;
function startLegacySensors() {
  let hadAny = false;
  if ('DeviceMotionEvent' in window) {
    window.addEventListener('devicemotion', (e) => {
      // acceleration (no gravity) preferred; else includingGravity
      const a = e.acceleration;
      const ag = e.accelerationIncludingGravity;
      if (a && a.x != null) {
        linAcc = { x: a.x || 0, y: a.y || 0, z: a.z || 0 };
      } else if (ag && ag.x != null) {
        accInc = { x: ag.x || 0, y: ag.y || 0, z: ag.z || 0 };
      }
      if (e.rotationRate) {
        const rr = e.rotationRate;
        // DeviceMotion rotationRate is in deg/s
        dmRotRate = { x: (rr.alpha || 0) * Math.PI/180,
                      y: (rr.beta || 0) * Math.PI/180,
                      z: (rr.gamma || 0) * Math.PI/180 };
      }
      dmTS = e.timeStamp;
    }, { passive: true });
    hadAny = true;
  }
  if ('DeviceOrientationEvent' in window) {
    window.addEventListener('deviceorientation', (e) => {
      // If available, build device->world quaternion from Euler angles
      if (e.alpha != null && e.beta != null && e.gamma != null) {
        // Build matrix then convert to quaternion (device->world)
        const a = (e.alpha||0) * Math.PI/180, b = (e.beta||0) * Math.PI/180, g = (e.gamma||0)*Math.PI/180;
        const cA=Math.cos(a), sA=Math.sin(a), cB=Math.cos(b), sB=Math.sin(b), cG=Math.cos(g), sG=Math.sin(g);
        // R = Rz(a) * Rx(b) * Ry(g)
        const m00 = cA*cG - sA*sB*sG;
        const m01 = -cB*sA;
        const m02 = cA*sG + cG*sA*sB;
        const m10 = cG*sA + cA*sB*sG;
        const m11 = cA*cB;
        const m12 = sA*sG - cA*cG*sB;
        const m20 = -cB*sG;
        const m21 = sB;
        const m22 = cB*cG;
        // Convert rotation matrix to quaternion
        const tr = m00 + m11 + m22;
        let qw, qx, qy, qz;
        if (tr > 0) {
          const S = Math.sqrt(tr + 1.0) * 2; // S=4*qw
          qw = 0.25 * S;
          qx = (m21 - m12) / S;
          qy = (m02 - m20) / S;
          qz = (m10 - m01) / S;
        } else if ((m00 > m11) && (m00 > m22)) {
          const S = Math.sqrt(1.0 + m00 - m11 - m22) * 2; // S=4*qx
          qw = (m21 - m12) / S;
          qx = 0.25 * S;
          qy = (m01 + m10) / S;
          qz = (m02 + m20) / S;
        } else if (m11 > m22) {
          const S = Math.sqrt(1.0 + m11 - m00 - m22) * 2; // S=4*qy
          qw = (m02 - m20) / S;
          qx = (m01 + m10) / S;
          qy = 0.25 * S;
          qz = (m12 + m21) / S;
        } else {
          const S = Math.sqrt(1.0 + m22 - m00 - m11) * 2; // S=4*qz
          qw = (m10 - m01) / S;
          qx = (m02 + m20) / S;
          qy = (m12 + m21) / S;
          qz = 0.25 * S;
        }
        qdw_abs = qNormalize([qw, qx, qy, qz]);
        usingAbsoluteOrientation = true;
        els.oriStatus.textContent = 'deviceorientation (fallback)';
      }
    }, true);
    hadAny = true;
  }
  return hadAny;
}

// Compute linear acceleration in device frame
function getLinearAccDev() {
  if (linAcc) return linAcc;
  if (accInc) {
    // Subtract gravity using quaternion
    // Need device->world quaternion
    let qdw;
    if (usingAbsoluteOrientation) {
      qdw = qdw_abs;
    } else {
      const qwd = Madgwick.getQwd();
      qdw = qConj(qwd);
    }
    // Gravity in world frame is [0,0,-G], rotate into device frame: g_dev = qwd ⊗ g_world ⊗ qwd*
    const qwd = usingAbsoluteOrientation ? qConj(qdw) : Madgwick.getQwd();
    const gDevVec = qRotateVec(qwd, { x: 0, y: 0, z: -G });
    return {
      x: (accInc.x || 0) - gDevVec.x,
      y: (accInc.y || 0) - gDevVec.y,
      z: (accInc.z || 0) - gDevVec.z,
    };
  }
  // Legacy DeviceMotion-only case: no data yet
  return { x: 0, y: 0, z: 0 };
}

// --------- Main loop ---------
function loop(t) {
  if (!running) return;
  if (lastFrameT == null) lastFrameT = t;
  const dt = Math.max(0, (t - lastFrameT) / 1000);
  lastFrameT = t;

  // Orientation fusion for legacy path using DeviceMotion gyro
  if (!usingAbsoluteOrientation && dmRotRate && dmTS) {
    const a = linAcc ?? accInc;
    if (a) {
      const dtDM = dt; // approximate
      Madgwick.update(dmRotRate.x, dmRotRate.y, dmRotRate.z, a.x, a.y, a.z, (mag?.x||0), (mag?.y||0), (mag?.z||0), dtDM);
      els.oriStatus.textContent = 'Madgwick (legacy)';
    }
  }

  // Get linear acceleration in device frame, apply small noise gate
  let aDev = getLinearAccDev();
  if (Math.abs(aDev.x) < GATE) aDev.x = 0;
  if (Math.abs(aDev.y) < GATE) aDev.y = 0;
  if (Math.abs(aDev.z) < GATE) aDev.z = 0;

  // Rotate to world frame using device->world quaternion
  let qdw;
  if (usingAbsoluteOrientation) qdw = qdw_abs;
  else qdw = qConj(Madgwick.getQwd());
  const aWorld = qRotateVec(qdw, aDev);

  // Integrate to world velocity
  vWorld = addScaled(vWorld, aWorld, dt);

  // Compute v magnitude and relativistic gamma with speed limit
  const vmag = norm3(vWorld.x, vWorld.y, vWorld.z);
  const c = getC();
  const vEff = limitSpeed(vmag, c);
  const gamma = gammaFromV(vEff, c);

  // Adaptive drift when near stillness (world frame magnitude)
  const amag = norm3(aWorld.x, aWorld.y, aWorld.z);
  if (amag === 0 || amag < GATE * 1.5) {
    const tau = getTau();
    const k = getK();
    const tauEff = tau / (1 + k * (gamma - 1));
    const decay = Math.exp(-dt / tauEff);
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

// --------- Controls & modal ----------
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
    lastFrameT = null;
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
