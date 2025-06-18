// The linear solver, and other functions for liquid simulation are taken from https://mikeash.com/pyblog/fluid-simulation-for-dummies.html
// with modifications to make it 2D and adapt to my usecase (f.e. renderExposure)


class FluidGrid {
  int N;          // Size of the grid (active cells)
  int iter;       // Solver iterations

  float dt;       // Timestep
  float diff;     // Diffusion rate
  float visc;     // Viscosity

  float[] s;      // Scratch space for density
  float[] density; // Density of "milk"

  float[] Vx;     // Velocity x-component
  float[] Vy;     // Velocity y-component
  float[] Vx0;    // Previous velocity x
  float[] Vy0;    // Previous velocity y

  FluidGrid(int N_val, float diffusion, float viscosity, float dt_val, int iter_val) {
    this.N = N_val;
    this.diff = diffusion;
    this.visc = viscosity;
    this.dt = dt_val;
    this.iter = iter_val;

    int gridSize = (N + 2) * (N + 2); // Total cells including boundaries
    this.s = new float[gridSize];
    this.density = new float[gridSize];

    this.Vx = new float[gridSize];
    this.Vy = new float[gridSize];
    this.Vx0 = new float[gridSize];
    this.Vy0 = new float[gridSize];
  }

  int IX(int x, int y) {
    x = constrain(x, 0, N + 1);
    y = constrain(y, 0, N + 1);
    return x + y * (N + 2);
  }

  void addDensity(int x, int y, float amount) {
    if (x >= 1 && x <= N && y >= 1 && y <= N) { // Only add to active grid cells
        this.density[IX(x, y)] += amount;
        this.density[IX(x, y)] = max(0, this.density[IX(x, y)]); // Optional: prevent negative density from accumulation if amount can be negative
    }
  }

  void setDensity(int x, int y, float value) {
    if (x >= 1 && x <= N && y >= 1 && y <= N) { // Only modify active grid cells
      this.density[IX(x, y)] = value;
    }
  }

  void addVelocity(int x, int y, float amountX, float amountY) {
     if (x >= 1 && x <= N && y >= 1 && y <= N) { // Only add to active grid cells
        int index = IX(x, y);
        this.Vx[index] += amountX;
        this.Vy[index] += amountY;
    }
  }

  void step() {
    // Store current velocities in Vx0, Vy0 to use as "previous" state for diffusion/advection
    System.arraycopy(Vx, 0, Vx0, 0, Vx.length);
    System.arraycopy(Vy, 0, Vy0, 0, Vy.length);
    System.arraycopy(density, 0, s, 0, density.length); // s becomes prev_density

    // Velocities
    diffuse(1, Vx0, Vx, visc, dt, iter, N); // Output to Vx0
    diffuse(2, Vy0, Vy, visc, dt, iter, N); // Output to Vy0

    project(Vx0, Vy0, Vx, Vy, iter, N); // Vx,Vy are p, div scratch. Output to Vx0, Vy0

    advect(1, Vx, Vx0, Vx0, Vy0, dt, N);  // Advect using projected Vx0,Vy0, output to Vx
    advect(2, Vy, Vy0, Vx0, Vy0, dt, N);  // Advect using projected Vx0,Vy0, output to Vy

    project(Vx, Vy, Vx0, Vy0, iter, N);   // Project current Vx,Vy. Output to Vx, Vy

    // Density
    // s already holds previous density state from System.arraycopy above
    diffuse(0, s, density, diff, dt, iter, N); // Diffuse current density, use s for scratch, output to s (new diffused state)
    advect(0, density, s, Vx, Vy, dt, N);     // Advect the diffused state (in s) using current velocities, output to density
    clipToCircle((N + 2)/2, (N + 2)/2, N/2);  
}

  // magic 
  void lin_solve(int b, float[] x, float[] x0, float a, float c, int iter_val, int N_active) {
    float cRecip = 1.0f / c;
    for (int k = 0; k < iter_val; k++) {
      for (int j = 1; j <= N_active; j++) {
        for (int i = 1; i <= N_active; i++) {
          x[IX(i, j)] =
            (x0[IX(i, j)]
            + a * (x[IX(i + 1, j)]
            + x[IX(i - 1, j)]
            + x[IX(i, j + 1)]
            + x[IX(i, j - 1)]
            )) * cRecip;
        }
      }
      set_bnd(b, x, N_active);
    }
  }

  void diffuse(int b, float[] x, float[] x0, float diff_rate, float dt_val, int iter_val, int N_active) {
    // x is the array to be filled (e.g. Vx0 from step), x0 is the source (e.g. Vx from step)
    float a = dt_val * diff_rate * (N_active * N_active);
    lin_solve(b, x, x0, a, 1 + 4 * a, iter_val, N_active);
  }

  void project(float[] velocX_in, float[] velocY_in, float[] p_scratch, float[] div_scratch, int iter_val, int N_active) {
    // velocX_in, velocY_in are the velocities to be projected (modified in place)
    // p_scratch, div_scratch are temporary arrays (passed as Vx, Vy from step())

    for (int j = 1; j <= N_active; j++) {
      for (int i = 1; i <= N_active; i++) {
        div_scratch[IX(i, j)] = -0.5f * (
          velocX_in[IX(i + 1, j)]
          - velocX_in[IX(i - 1, j)]
          + velocY_in[IX(i, j + 1)]
          - velocY_in[IX(i, j - 1)]
          ) / N_active;
        p_scratch[IX(i, j)] = 0;
      }
    }
    set_bnd(0, div_scratch, N_active);
    set_bnd(0, p_scratch, N_active);
    lin_solve(0, p_scratch, div_scratch, 1, 4, iter_val, N_active);

    for (int j = 1; j <= N_active; j++) {
      for (int i = 1; i <= N_active; i++) {
        velocX_in[IX(i, j)] -= 0.5f * (p_scratch[IX(i + 1, j)] - p_scratch[IX(i - 1, j)]) * N_active;
        velocY_in[IX(i, j)] -= 0.5f * (p_scratch[IX(i, j + 1)] - p_scratch[IX(i, j - 1)]) * N_active;
      }
    }
    set_bnd(1, velocX_in, N_active);
    set_bnd(2, velocY_in, N_active);
  }

  void advect(int b, float[] d_out, float[] d_in, float[] velocX_field, float[] velocY_field, float dt_val, int N_active) {
    // d_out is the array to be filled (e.g. Vx from step), d_in is the source (e.g. Vx0 from step)
    float i0_f, i1_f, j0_f, j1_f; // _f for float
    float dtx = dt_val * N_active;
    float dty = dt_val * N_active;

    float s0, s1, t0, t1;
    float tmp1, tmp2, current_x_coord_f, current_y_coord_f; // _f for float

    float Nfloat = N_active;
    float ifloat, jfloat; // Loop counters as floats
    // int i, j; // Loop counters as ints (already parameter or local)

    for (int j_int = 1, j = 1; j_int <= N_active; j_int++, j++) {
      jfloat = j;
      for (int i_int = 1, i = 1; i_int <= N_active; i_int++, i++) {
        ifloat = i;
        tmp1 = dtx * velocX_field[IX(i_int, j_int)];
        tmp2 = dty * velocY_field[IX(i_int, j_int)];
        current_x_coord_f = ifloat - tmp1;
        current_y_coord_f = jfloat - tmp2;

        if (current_x_coord_f < 0.5f) current_x_coord_f = 0.5f;
        if (current_x_coord_f > Nfloat + 0.5f) current_x_coord_f = Nfloat + 0.5f;
        i0_f = floor(current_x_coord_f);
        i1_f = i0_f + 1.0f;

        if (current_y_coord_f < 0.5f) current_y_coord_f = 0.5f;
        if (current_y_coord_f > Nfloat + 0.5f) current_y_coord_f = Nfloat + 0.5f;
        j0_f = floor(current_y_coord_f);
        j1_f = j0_f + 1.0f;

        s1 = current_x_coord_f - i0_f;
        s0 = 1.0f - s1;
        t1 = current_y_coord_f - j0_f;
        t0 = 1.0f - t1;

        int i0i = (int)i0_f;
        int i1i = (int)i1_f;
        int j0i = (int)j0_f;
        int j1i = (int)j1_f;

        d_out[IX(i_int, j_int)] =
          s0 * (t0 * d_in[IX(i0i, j0i)] + t1 * d_in[IX(i0i, j1i)]) +
          s1 * (t0 * d_in[IX(i1i, j0i)] + t1 * d_in[IX(i1i, j1i)]);
      }
    }
    set_bnd(b, d_out, N_active);
  }

  void set_bnd(int b, float[] x, int N_active) {
    for (int i = 1; i <= N_active; i++) {
      x[IX(0, i)]     = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
      x[IX(N_active + 1, i)] = b == 1 ? -x[IX(N_active, i)] : x[IX(N_active, i)];
      x[IX(i, 0)]     = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
      x[IX(i, N_active + 1)] = b == 2 ? -x[IX(i, N_active)] : x[IX(i, N_active)];
    }

    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N_active + 1)] = 0.5f * (x[IX(1, N_active + 1)] + x[IX(0, N_active)]);
    x[IX(N_active + 1, 0)] = 0.5f * (x[IX(N_active, 0)] + x[IX(N_active + 1, 1)]);
    x[IX(N_active + 1, N_active + 1)] = 0.5f * (x[IX(N_active, N_active + 1)] + x[IX(N_active + 1, N_active)]);
  }

  void renderD(PGraphics pg, int r_scale, float exposure) {
    pg.noStroke();
    for (int j = 1; j <= N; j++) {
      for (int i = 1; i <= N; i++) {
        float d = density[IX(i, j)];
        if (d > 0.01f) {
          float lerp = constrain(d / exposure, 0, 1);
          pg.fill(lerpColor(coffeeColor, milkColor, lerp));
          pg.rect(i * r_scale, j * r_scale, r_scale, r_scale);
        }
      }
    }
  }

  void reset() {
    int size = (N + 2) * (N + 2);
    for (int i = 0; i < size; i++) {
      density[i] = s[i] = 0;
      Vx[i] = Vy[i] = Vx0[i] = Vy0[i] = 0;
    }
  }


  void clipToCircle(int cx, int cy, int rad) {
    int r2 = rad * rad;
    for (int j = 1; j <= N; j++) {
      for (int i = 1; i <= N; i++) {
        int d2 = (i - cx)*(i - cx) + (j - cy)*(j - cy);
        if (d2 > r2) {                // outside mug
          int idx = IX(i,j);
          density[idx] = 0;
          Vx[idx] = Vy[idx] = 0;
        }
      }
    }
  }
}
