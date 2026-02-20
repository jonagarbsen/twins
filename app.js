// Twin Timer with world-frame integration
// 1) Use DeviceMotion for acceleration (prefer gravity-free).
// 2) Use DeviceOrientation (alpha, beta, gamma) to build rotation R_device->world.
// 3) Rotate accel to world frame, integrate to world velocity.
// 4) Drift-to-zero when nearly still; limit speed to < c; update proper time.

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
};

let running = false;
let lastFrameT = null;

// Acceleration (device frame, m/s^2)
let aDev = { x: 0, y: 0, z: 0 };
// Velocity (world frame, m/s)
let vWorld = { x: 0, y: 0, z: 0 };

// Orientation matrix R (device -> world)
let R = identity3();

// High-pass filter memory for gravity removal fallback
const hp = { x: 0, y: 0, z: 0, lastX: null, lastY: null, lastZ: null };
const HP_ALPHA = 0.98;

// Proper time (s)
let tauVal = 0;

// Settings
const GATE = 0.05; // quiet gate for tiny accel

// -------------- Math helpers --------------
function identity3() {
  return [1,0,0, 0,1,0, 0,0,1];
}
function mult3x3Vec3(M, v) {
  return {
    x: M[0]*v.x + M[1]*v.y + M[2]*v.z,
    y: M[3]*v.x + M[4]*v.y + M[5]*v.z,
    z: M[6]*v.x + M[7]*v.y + M[8]*v.z,
  };
}
function norm3(x, y, z) { return Math.hypot(x, y, z); }
function addScaled(v, a, dt) { return { x: v.x + a.x * dt, y: v.y + a.y * dt, z: v.z + a.z * dt }; }
function formatTau(sec) {
  const ms = Math.floor((sec % 1) * 1000);
  const total = Math.floor(sec);
  const s = total % 60;
  const m = Math.floor((total / 60) % 60);
  const h = Math.floor(total / 3600);
  const pad = (n, w = 2) => n.toString().padStart(w, '0');
  return `${pad(h)}:${pad(m)}:${pad(s)}.${pad(ms, 3)}`;
}
function getC() {
  const c = parseFloat(els.cval.value);
  return Number.isFinite(c) && c > 0 ? c : 1.0;
}
function getTau() {
  const t = parseFloat(els.tau.value);
  return Number.isFinite(t) && t > 0 ? t : 3.0;
}
function limitSpeed(vmag, c) {
  return c * Math.tanh(vmag / c);
}
function gammaFromV(vmag, c) {
  const beta2 = (vmag * vmag) / (c * c);
  if (beta2 <= 0) return 1;
  if (beta2 >= 1) return 1e12;
  return 1 / Math.sqrt(1 - beta2);
}

// -------------- Orientation handling --------------
// Build R from DeviceOrientation angles (alpha, beta, gamma), Z–X′–Y″ convention.
// alpha: z (compass), beta: x (pitch), gamma: y (roll).
function eulerToMatrix(alphaRad, betaRad, gammaRad) {
  const cA = Math.cos(alphaRad), sA = Math.sin(alphaRad);
  const cB = Math.cos(betaRad), sB = Math.sin(betaRad);
  const cG = Math.cos(gammaRad), sG = Math.sin(gammaRad);

  // R = Rz(alpha) * Rx(beta) * Ry(gamma)
  const m00 = cA*cG - sA*sB*sG;
  const m01 = -cB*sA;
  const m02 = cA*sG + cG*sA*sB;

  const m10 = cG*sA + cA*sB*sG;
  const m11 = cA*cB;
  const m12 = sA*sG - cA*cG*sB;

  const m20 = -cB*sG;
  const m21 = sB;
  const m22 = cB*cG;

  return [m00,m01,m02, m10,m11,m12, m20,m21,m22];
}

// Adjust for screen orientation (portrait vs landscape).
// For simplicity, we assume portrait-primary. If landscape, rotate around z by ±90°.
function applyScreenOrientation(M) {
  const type = (screen.orientation && screen.orientation.type) || '';
  if (type.startsWith('landscape')) {
    const angle = (screen.orientation.angle || (type.endsWith('-secondary') ? 270 : 90)) * Math.PI / 180;
    const c = Math.cos(angle), s = Math.sin(angle);
    const Rz = [c,-s,0, s,c,0, 0,0,1];
    return mult3x3(Rz, M);
  }
  return M;
}
function mult3x3(A, B) {
  const r = new Array(9);
  r[0]=A[0]*B[0]+A[1]*B[3]+A[2]*B[6];
  r[1]=A[0]*B[1]+A[1]*B[4]+A[2]*B[7];
  r[2]=A[0]*B[2]+A[1]*B[5]+A[2]*B[8];
  r[3]=A[3]*B[0]+A[4]*B[3]+A[5]*B[6];
  r[4]=A[3]*B[1]+A[4]*B[4]+A[5]*B[7];
  r[5]=A[3]*B[2]+A[4]*B[5]+A[5]*B[8];
  r[6]=A[6]*B[0]+A[7]*B[3]+A[8]*B[6];
  r[7]=A[6]*B[1]+A[7]*B[4]+A[8]*B[7];
  r[8]=A[6]*B[2]+A[7]*B[5]+A[8]*B[8];
  return r;
}

function startOrientation() {
  if (!('DeviceOrientationEvent' in window)) {
    els.oriStatus.textContent = 'unsupported';
    return false;
  }
  // Some platforms expose deviceorientationabsolute; use it if available.
  const handler = (e) => {
    // alpha (0..360), beta (-180..180), gamma (-90..90)
    let a = e.alpha, b = e.beta, g = e.gamma;
    if (a == null || b == null || g == null) {
      els.oriStatus.textContent = 'no data';
      return;
    }
    // Convert to radians
    const aRad = (a || 0) * Math.PI / 180;
    const bRad = (b || 0) * Math.PI / 180;
    const gRad = (g || 0) * Math.PI / 180;

    let M = eulerToMatrix(aRad, bRad, gRad);
    M = applyScreenOrientation(M);
    R = M;
    els.oriStatus.textContent = 'ok';
  };

  // Prefer absolute orientation if available
  window.addEventListener('deviceorientationabsolute', handler, true);
  window.addEventListener('deviceorientation', handler, true);
  return true;
}

// -------------- Sensors (acceleration) --------------
function startAccel() {
  if (!('DeviceMotionEvent' in window)) {
    return false;
  }
  window.addEventListener('devicemotion', onMotion, { passive: true });
  return true;
}

function onMotion(e) {
  // Prefer acceleration without gravity
  let ax = e.acceleration?.x ?? null;
  let ay = e.acceleration?.y ?? null;
  let az = e.acceleration?.z ?? null;

  if (ax == null || ay == null || az == null) {
    // Fallback: accelerationIncludingGravity with high-pass
    const gx = e.accelerationIncludingGravity?.x ?? 0;
    const gy = e.accelerationIncludingGravity?.y ?? 0;
    const gz = e.accelerationIncludingGravity?.z ?? 0;

    hp.x = HP_ALPHA * (hp.x + gx - (hp.lastX ?? gx));
    hp.y = HP_ALPHA * (hp.y + gy - (hp.lastY ?? gy));
    hp.z = HP_ALPHA * (hp.z + gz - (hp.lastZ ?? gz));
    hp.lastX = gx; hp.lastY = gy; hp.lastZ = gz;

    ax = hp.x; ay = hp.y; az = hp.z;
  }

  // Noise gate in device frame
  if (Math.abs(ax) < GATE) ax = 0;
  if (Math.abs(ay) < GATE) ay = 0;
  if (Math.abs(az) < GATE) az = 0;

  aDev = { x: ax, y: ay, z: az };
}

// -------------- Main loop --------------
function loop(t) {
  if (!running) return;
  if (lastFrameT == null) lastFrameT = t;
  const dt = Math.max(0, (t - lastFrameT) / 1000);
  lastFrameT = t;

  // Rotate device accel to world frame
  const aWorld = mult3x3Vec3(R, aDev);

  // Integrate acceleration to world velocity
  vWorld = addScaled(vWorld, aWorld, dt);

  // Drift-to-zero when nearly still (use world-frame magnitude)
  const amag = norm3(aWorld.x, aWorld.y, aWorld.z);
  if (amag === 0 || amag < GATE * 1.5) {
    const decay = Math.exp(-dt / getTau());
    vWorld.x *= decay; vWorld.y *= decay; vWorld.z *= decay;
  }

  const vmag = norm3(vWorld.x, vWorld.y, vWorld.z);
  const c = getC();
  const vEff = limitSpeed(vmag, c);
  const gamma = gammaFromV(vEff, c);

  // Proper time update: Δτ = Δt / γ
  tauVal += dt / gamma;

  // UI updates (show world-frame accel magnitude and speed)
  els.timer.textContent = formatTau(tauVal);
  els.accel.textContent = `${aWorld.x.toFixed(3)}, ${aWorld.y.toFixed(3)}, ${aWorld.z.toFixed(3)}`;
  els.vraw.textContent = vmag.toFixed(3);
  els.veff.textContent = vEff.toFixed(3);
  els.gamma.textContent = gamma.toFixed(6);

  requestAnimationFrame(loop);
}

// -------------- Controls & modal --------------
function openSettings() {
  els.modalBackdrop.hidden = false;
  if (typeof els.settings.showModal === 'function') {
    els.settings.showModal();
  } else {
    els.settings.style.display = 'block';
  }
}
function closeSettings() {
  els.modalBackdrop.hidden = true;
  if (typeof els.settings.close === 'function') {
    els.settings.close();
  } else {
    els.settings.style.display = 'none';
  }
}

els.openSettings.addEventListener('click', openSettings);
els.closeSettings.addEventListener('click', closeSettings);
els.closeSettingsFooter.addEventListener('click', closeSettings);
els.modalBackdrop.addEventListener('click', closeSettings);

els.startStop.addEventListener('click', () => {
  if (!running) {
    startOrientation();
    startAccel();
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
  aDev = { x: 0, y: 0, z: 0 };
  hp.x = hp.y = hp.z = 0; hp.lastX = hp.lastY = hp.lastZ = null;
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
