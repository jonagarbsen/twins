// Twin Paradox Timer
// Integrates acceleration (no gravity if available) to velocity,
// applies drift-to-zero for stillness, limits speed to < c,
// and updates proper time tau = ∫ dt / gamma.

const els = {
  timer: document.getElementById('timer'),
  startStop: document.getElementById('startStop'),
  reset: document.getElementById('reset'),
  status: document.getElementById('status'),
  accel: document.getElementById('accel'),
  vraw: document.getElementById('vraw'),
  veff: document.getElementById('veff'),
  gamma: document.getElementById('gamma'),
  tau: document.getElementById('tau'),
  gate: document.getElementById('gate'),
  cval: document.getElementById('cval'),
  limiter: document.getElementById('limiter'),
};

let running = false;
let lastFrameT = null;

// Latest acceleration sample (m/s^2), device coordinates
let a = { x: 0, y: 0, z: 0 };
// Current velocity estimate (m/s), device coordinates
let v = { x: 0, y: 0, z: 0 };

// High-pass filter memory for gravity removal fallback
const hp = { x: 0, y: 0, z: 0, lastX: null, lastY: null, lastZ: null };
const HP_ALPHA = 0.98;

// Proper time accumulator (seconds)
let tau = 0;

// Helper math
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
function getTau() {
  const t = parseFloat(els.tau.value);
  return Number.isFinite(t) && t > 0 ? t : 3.0;
}
function getGate() {
  const g = parseFloat(els.gate.value);
  return Number.isFinite(g) && g >= 0 ? g : 0.05;
}
function getC() {
  const c = parseFloat(els.cval.value);
  return Number.isFinite(c) && c > 0 ? c : 1.0;
}
function limitSpeed(vmag, c, mode) {
  if (vmag <= 0) return 0;
  if (mode === 'tanh') {
    // v_eff = c * tanh(v/c) -> strictly less than c
    return c * Math.tanh(vmag / c);
  }
  // clip to 0.999c
  return Math.min(vmag, 0.999 * c);
}
function gammaFromV(vmag, c) {
  const beta2 = (vmag * vmag) / (c * c);
  if (beta2 >= 1) return 1e12; // safeguard
  return 1 / Math.sqrt(1 - beta2);
}

// Motion sensor setup
function startSensors() {
  if (!('DeviceMotionEvent' in window)) {
    els.status.textContent = 'DeviceMotion not supported on this browser.';
    return false;
  }
  // Android Chrome usually needs only a user gesture; no explicit permission.
  window.addEventListener('devicemotion', onMotion, { passive: true });
  els.status.textContent = 'listening to devicemotion…';
  return true;
}

function onMotion(e) {
  // Prefer acceleration (without gravity)
  let ax = e.acceleration && e.acceleration.x != null ? e.acceleration.x : null;
  let ay = e.acceleration && e.acceleration.y != null ? e.acceleration.y : null;
  let az = e.acceleration && e.acceleration.z != null ? e.acceleration.z : null;

  if (ax == null || ay == null || az == null) {
    // Fall back to includingGravity with a high-pass filter
    const gx = e.accelerationIncludingGravity?.x ?? 0;
    const gy = e.accelerationIncludingGravity?.y ?? 0;
    const gz = e.accelerationIncludingGravity?.z ?? 0;

    hp.x = HP_ALPHA * (hp.x + gx - (hp.lastX ?? gx));
    hp.y = HP_ALPHA * (hp.y + gy - (hp.lastY ?? gy));
    hp.z = HP_ALPHA * (hp.z + gz - (hp.lastZ ?? gz));
    hp.lastX = gx; hp.lastY = gy; hp.lastZ = gz;

    ax = hp.x; ay = hp.y; az = hp.z;
    els.status.textContent = 'using high-pass filter (no gravity-free accel)';
  } else {
    els.status.textContent = 'using acceleration (no gravity)';
  }

  // Simple noise gate
  const G = getGate();
  if (Math.abs(ax) < G) ax = 0;
  if (Math.abs(ay) < G) ay = 0;
  if (Math.abs(az) < G) az = 0;

  a = { x: ax, y: ay, z: az };
  els.accel.textContent = `${ax.toFixed(3)}, ${ay.toFixed(3)}, ${az.toFixed(3)}`;
}

// Main loop
function loop(t) {
  if (!running) return;
  if (lastFrameT == null) lastFrameT = t;
  const dt = Math.max(0, (t - lastFrameT) / 1000); // seconds
  lastFrameT = t;

  // Integrate acceleration to velocity
  v = addScaled(v, a, dt);

  // Drift-to-zero when near stillness
  const amag = norm3(a.x, a.y, a.z);
  const tauDecay = getTau();
  if (amag === 0 || amag < getGate() * 1.5) {
    // Exponential decay: v <- v * exp(-dt/τ)
    const decay = Math.exp(-dt / tauDecay);
    v.x *= decay; v.y *= decay; v.z *= decay;
  }

  const vmag = norm3(v.x, v.y, v.z);
  const c = getC();
  const vEff = limitSpeed(vmag, c, els.limiter.value);
  const gamma = gammaFromV(vEff, c);

  // Proper time update: dτ = dt / γ
  const dTau = dt / gamma;
  tau += dTau;

  // Update UI
  els.timer.textContent = formatTau(tau);
  els.vraw.textContent = vmag.toFixed(3);
  els.veff.textContent = vEff.toFixed(3);
  els.gamma.textContent = gamma.toFixed(6);

  requestAnimationFrame(loop);
}

// Controls
els.startStop.addEventListener('click', async () => {
  if (!running) {
    // Start sensors on user gesture
    const ok = startSensors();
    if (!ok) return;
    running = true;
    lastFrameT = null;
    requestAnimationFrame(loop);
    els.startStop.textContent = 'Pause';
  } else {
    running = false;
    els.startStop.textContent = 'Start';
    els.status.textContent = 'paused';
  }
});

els.reset.addEventListener('click', () => {
  // Reset timer and speed
  tau = 0;
  v = { x: 0, y: 0, z: 0 };
  a = { x: 0, y: 0, z: 0 };
  hp.x = hp.y = hp.z = 0; hp.lastX = hp.lastY = hp.lastZ = null;
  els.timer.textContent = formatTau(0);
  els.vraw.textContent = '0.00';
  els.veff.textContent = '0.00';
  els.gamma.textContent = '1.000000';
  els.status.textContent = 'reset';
});

document.addEventListener('visibilitychange', () => {
  // Pause when backgrounded to avoid weird timing
  if (document.hidden && running) {
    running = false;
    els.startStop.textContent = 'Start';
    els.status.textContent = 'paused (backgrounded)';
  }
});