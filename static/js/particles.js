/**
 * Ambient Particle Background — Premium floating particles with connections
 * Creates a subtle, interactive starfield/DNA-helix inspired animation
 */
(function () {
    const canvas = document.getElementById('particle-canvas');
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    let W, H;
    let particles = [];
    let animFrame;
    const PARTICLE_COUNT = 55;
    const CONNECTION_DIST = 140;
    const MOUSE = { x: -1000, y: -1000 };

    function resize() {
        W = window.innerWidth;
        H = window.innerHeight;
        canvas.width = W;
        canvas.height = H;
    }

    function createParticle() {
        const colors = [
            'rgba(255, 146, 51, 0.35)',  // orange
            'rgba(255, 192, 69, 0.30)',  // gold
            'rgba(99, 102, 241, 0.25)',  // indigo
            'rgba(0, 212, 255, 0.20)',   // cyan
            'rgba(78, 159, 61, 0.20)',   // green
        ];
        return {
            x: Math.random() * W,
            y: Math.random() * H,
            vx: (Math.random() - 0.5) * 0.3,
            vy: (Math.random() - 0.5) * 0.3,
            r: Math.random() * 2 + 0.5,
            color: colors[Math.floor(Math.random() * colors.length)],
            pulse: Math.random() * Math.PI * 2,
            pulseSpeed: 0.005 + Math.random() * 0.015,
        };
    }

    function init() {
        resize();
        particles = [];
        for (let i = 0; i < PARTICLE_COUNT; i++) {
            particles.push(createParticle());
        }
    }

    function draw() {
        ctx.clearRect(0, 0, W, H);

        // Update & draw particles
        particles.forEach(p => {
            p.x += p.vx;
            p.y += p.vy;
            p.pulse += p.pulseSpeed;

            // Wrap around edges
            if (p.x < -10) p.x = W + 10;
            if (p.x > W + 10) p.x = -10;
            if (p.y < -10) p.y = H + 10;
            if (p.y > H + 10) p.y = -10;

            // Mouse repulsion (subtle)
            const dx = p.x - MOUSE.x;
            const dy = p.y - MOUSE.y;
            const dist = Math.sqrt(dx * dx + dy * dy);
            if (dist < 120) {
                p.vx += (dx / dist) * 0.02;
                p.vy += (dy / dist) * 0.02;
            }

            // Dampen velocity
            p.vx *= 0.999;
            p.vy *= 0.999;

            // Pulsing size
            const r = p.r + Math.sin(p.pulse) * 0.4;

            ctx.beginPath();
            ctx.arc(p.x, p.y, Math.max(r, 0.3), 0, Math.PI * 2);
            ctx.fillStyle = p.color;
            ctx.fill();
        });

        // Draw connections
        for (let i = 0; i < particles.length; i++) {
            for (let j = i + 1; j < particles.length; j++) {
                const dx = particles[i].x - particles[j].x;
                const dy = particles[i].y - particles[j].y;
                const dist = Math.sqrt(dx * dx + dy * dy);
                if (dist < CONNECTION_DIST) {
                    const alpha = (1 - dist / CONNECTION_DIST) * 0.08;
                    ctx.strokeStyle = `rgba(255, 146, 51, ${alpha})`;
                    ctx.lineWidth = 0.5;
                    ctx.beginPath();
                    ctx.moveTo(particles[i].x, particles[i].y);
                    ctx.lineTo(particles[j].x, particles[j].y);
                    ctx.stroke();
                }
            }
        }

        animFrame = requestAnimationFrame(draw);
    }

    // Mouse tracking
    document.addEventListener('mousemove', (e) => {
        MOUSE.x = e.clientX;
        MOUSE.y = e.clientY;
    });

    document.addEventListener('mouseleave', () => {
        MOUSE.x = -1000;
        MOUSE.y = -1000;
    });

    window.addEventListener('resize', () => {
        resize();
    });

    init();
    draw();
})();
