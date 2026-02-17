// Twin Timer (clean UI)
// Proper time tau updates via SR time dilation using accelerometer-derived speed.
// v_eff = c * tanh(v_raw / c) ensures speed < c.

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

  cval: document.getElementById('cval'),
  tau: document.getElementById('tau'),
};

let running = false;
let lastFrameT = null;

// Acceleration (m/s^2)
let a = { x: 0, y: 0, z: 0 };
// Velocity (m/s)
let v = { x: 0, y: 0, z: 0 };

// High-pass filter memory for gravity removal fallback
const hp = { x: 0, y: 0, z: 0, lastX: null, lastY: null, lastZ: null };
const HP_ALPHA = 0.98;

// Proper time (s)
let tauVal = 0;

// Constants / settings
const GATE = 0.05; // quiet gate for tiny accel
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

// Config getters
function getC() {
  const c = parseFloat(els.cval.value);
  return Number.isFinite(c) && c > 0 ? c : 1.0;
}
function getTau() {
  const t = parseFloat(els.tau.value);
  return Number.isFinite(t) && t > 0 ? t : 3.0;
}

// Speed limiting
function limitSpeed(vmag, c) {
  return c * Math.tanh(vmag / c); // strictly less than c
}
function gammaFromV(vmag, c) {
  const beta2 = (vmag * vmag) / (c * c);
  if (beta2 <= 0) return 1;
  if (beta2 >= 1) return 1e12; // guard
  return 1 / Math.sqrt(1 - beta2);
}

// Sensors
function startSensors() {
  if (!('DeviceMotionEvent' in window)) {
    // Minimal UI: silently ignore on unsupported browsers
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

  // Noise gate
  if (Math.abs(ax) < GATE) ax = 0;
  if (Math.abs(ay) < GATE) ay = 0;
  if (Math.abs(az) < GATE) az = 0;

  a = { x: ax, y: ay, z: az };
  // Update settings panel live
  els.accel.textContent = `${ax.toFixed(3)}, ${ay.toFixed(3)}, ${az.toFixed(3)}`;
}

// Main loop
function loop(t) {
  if (!running) return;
  if (lastFrameT == null) lastFrameT = t;
  const dt = Math.max(0, (t - lastFrameT) / 1000);
  lastFrameT = t;

  // Integrate acceleration to velocity
  v = addScaled(v, a, dt);

  // Drift-to-zero when nearly still
  const amag = norm3(a.x, a.y, a.z);
  if (amag === 0 || amag < GATE * 1.5) {
    const decay = Math.exp(-dt / getTau());
    v.x *= decay; v.y *= decay; v.z *= decay;
  }

  const vmag = norm3(v.x, v.y, v.z);
  const c = getC();
  const vEff = limitSpeed(vmag, c);
  const gamma = gammaFromV(vEff, c);

  // Proper time update
  tauVal += dt / gamma;

  // UI
  els.timer.textContent = formatTau(tauVal);
  els.vraw.textContent = vmag.toFixed(3);
  els.veff.textContent = vEff.toFixed(3);
  els.gamma.textContent = gamma.toFixed(6);

  requestAnimationFrame(loop);
}

// Controls
function openSettings() {
  els.modalBackdrop.hidden = false;
  if (typeof els.settings.showModal === 'function') {
    els.settings.showModal();
  } else {
    // Fallback for browsers without <dialog> support
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
    if (!startSensors()) {
      // If sensors are unavailable, still run the timer (non-relativistic)
    }
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
  v = { x: 0, y: 0, z: 0 };
  a = { x: 0, y: 0, z: 0 };
  hp.x = hp.y = hp.z = 0; hp.lastX = hp.lastY = hp.lastZ = null;
  els.timer.textContent = formatTau(0);
  els.vraw.textContent = '0.00';
  els.veff.textContent = '0.00';
  els.gamma.textContent = '1.000000';
  // keep running state unchanged
});

document.addEventListener('visibilitychange', () => {
  if (document.hidden && running) {
    running = false;
    els.startStop.textContent = 'Start';
  }
});
