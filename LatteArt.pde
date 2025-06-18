import controlP5.*;

PGraphics drawLayer;
PGraphics UILayer;
String exportFolder = sketchPath("exports");
ControlP5 cp5;
PImage mugOverlay;
FluidGrid fluid;
boolean overlayLoaded = false;
int simOffsetX = 66;
int simOffsetY = 137;

// Simulation parameters
int    N = 256;      // grid resolution
float  dt = 0.28f;    // time-step
int    ITER = 20;
float  currentMilkAmount;
float  currentForceAmount;
float  currentViscosity;
float  currentDiffusion;
int    brushSize;
float  renderExposure;
int    SCALE = 1; // 1 px per cell

color coffeeColor = color(100, 65, 23);
color milkColor   = color(250, 245, 230);

int mugRadius   = N / 2;
int mugCX       = (N + 2) / 2;
int mugCY       = mugCX;

Slider milkAmountSlider,
  forceAmountSlider,
  viscositySlider,
  diffusionSlider,
  exposureSlider,
  brushSizeSlider;

// buttons
boolean paused = false;      // current state
Toggle  pauseToggle;
boolean symmetryEnabled = false;
Toggle symmetryToggle;
boolean showAbout = false;
Button  closeAboutBtn;

void settings() {
  size(700, 512);
}

void setup() {
  currentMilkAmount  = 70;
  currentForceAmount = 0.0015f;
  currentViscosity   = 1f;
  currentDiffusion   = 0.0f;
  brushSize          = 8;
  renderExposure     = 95;

  mugOverlay = loadImage("data/overlay.png");
  overlayLoaded = (mugOverlay != null);
  fluid = new FluidGrid(N,
    currentDiffusion,
    currentViscosity,
    dt,
    ITER);

  drawLayer = createGraphics(700, 512);
  UILayer = createGraphics(700, 512);
  cp5 = new ControlP5(this);

  int uiX = (N + 2) * SCALE + 190;
  int uiY = 50;
  int w   = 160;
  int h   = 20;
  int step = 30;

  milkAmountSlider = cp5.addSlider("currentMilkAmount")
    .setPosition(uiX, uiY)
    .setSize(w, h)
    .setRange(1, 450)
    .setValue(currentMilkAmount)
    .setLabel("Milk Amount");

  uiY += step;
  forceAmountSlider = cp5.addSlider("currentForceAmount")
    .setPosition(uiX, uiY)
    .setSize(w, h)
    .setRange(0f, 0.03f)
    .setDecimalPrecision(4)
    .setValue(currentForceAmount)
    .setLabel("Force Amount");


  uiY += step;
  viscositySlider = cp5.addSlider("currentViscosity")
    .setPosition(uiX, uiY)
    .setSize(w, h)
    .setRange(0f, 1f)
    .setDecimalPrecision(4)
    .setValue(currentViscosity)
    .setLabel("Viscosity");

  uiY += step;
  diffusionSlider = cp5.addSlider("currentDiffusion")
    .setPosition(uiX, uiY)
    .setSize(w, h)
    .setRange(0.0f, 0.0001f)
    .setDecimalPrecision(6)
    .setValue(currentDiffusion)
    .setLabel("Milk Diffusion");

  uiY += step;
  exposureSlider  = cp5.addSlider("renderExposure")
    .setPosition(uiX, uiY)
    .setSize(w, h)
    .setRange(0, 98)
    .setValue(renderExposure)
    .setLabel("Render Exposure");

  uiY += step;
  brushSizeSlider  = cp5.addSlider("brushSize")
    .setPosition(uiX, uiY)
    .setSize(w, h)
    .setRange(5, 20)
    .setValue(brushSize)
    .setLabel("Brush Size");

  uiY += step * 2;
  cp5.addButton("resetFluidButton")
    .setPosition(uiX, uiY)
    .setSize(w, h + 10)
    .setLabel("Reset Milk")
    .setColorLabel(color(255))
    .setColorBackground(color(223, 83, 61))
    .setColorForeground(color(234, 137, 122));

  uiY += step * 2;
  pauseToggle = cp5.addToggle("paused")
    .setPosition(uiX, uiY)
    .setSize(w/4, h)
    .setMode(ControlP5.SWITCH)
    .setLabel("Pause")
    .setColorBackground(color(1, 173, 254))
    .setColorForeground(color(0, 46, 92))
    .setColorActive(color(0, 46, 92));

  uiX += 80;
  symmetryToggle = cp5.addToggle("symmetryEnabled")
    .setPosition(uiX, uiY)
    .setSize(w/4, h)
    .setValue(false)
    .setMode(ControlP5.SWITCH)
    .setLabel("Symmetry")
    .setColorBackground(color(1, 173, 254))
    .setColorForeground(color(0, 46, 92))
    .setColorActive(color(0, 46, 92));

  uiX -= 80;
  uiY += step * 2;
  cp5.addButton("exportPNGButton")
    .setPosition(uiX, uiY)
    .setSize(w, h + 10)
    .setLabel("Export PNG")
    .setColorForeground(color(111, 202, 108))
    .setColorBackground(color(91, 219, 87))
    .setColorLabel(color(0));

  uiY += step * 2;
  cp5.addButton("aboutButton")
    .setPosition(uiX, uiY)
    .setSize(w, h + 10)
    .setLabel("About");

  PFont myFont = createFont("Control Font", 12);
  Group aboutGroup = cp5.addGroup("aboutGroup")
    .setPosition(0, 0)
    .setSize(width, height)
    .setBackgroundColor(color(0, 160))
    .setBackgroundHeight(height)
    .hide();

  cp5.addTextarea("aboutText")
    .setGroup(aboutGroup)
    .setPosition(width/2 - 160, height/2 - 110)
    .setSize(320, 180)
    .setText("Latte Art Simulator\n"
    + "----------------------------------\n"
    + "Left mouse  –  pour milk\n"
    + "Right mouse –  erase milk\n"
    + "Pause btn   –  freezes/unfreezes simulation\n"
    + "Symmetry    -  toggles symmetry aroud Y axis\n"
    + "Reset milk  -  clears milk and resets flow grid\n"
    + "Sliders     -  as labels say\n"
    + "NOTE: if you increase brush size, it's recommended to decrease force and icrease milk amount\n"
    + "Export      –  exports PNG to (exports/)\n\n"
    + "© 2025 Timotej Dzugas")
    .setFont(myFont)
    .setColorBackground(color(0, 0))
    .setColorForeground(color(0, 0))
    .setColor(color(255));

  cp5.addButton("closeAboutBtn")
    .setGroup(aboutGroup)
    .setPosition(width/2 - 40, height/2 + 90)
    .setSize(80, 25)
    .setLabel("Close");
  println("Setup done.");
}

public void resetFluidButton() {
  fluid.reset();
  println("Fluid simulation reset.");
}

public void handleBrush() {
  boolean overSim =
    mouseX >= simOffsetX &&
    mouseX <  simOffsetX + drawLayer.width  &&
    mouseY >= simOffsetY &&
    mouseY <  simOffsetY + drawLayer.height;

  if (overSim && mousePressed) {
    int baseX = floor((mouseX - simOffsetX) / (float)SCALE);
    int baseY = floor((mouseY - simOffsetY) / (float)SCALE);
    int r2    = mugRadius * mugRadius;

    for (int dy = -brushSize; dy <= brushSize; dy++) {
      for (int dx = -brushSize; dx <= brushSize; dx++) {

        // skip offsets outside the circular mug
        int gx = baseX + dx - mugCX;
        int gy = baseY + dy - mugCY;
        if (gx*gx + gy*gy > r2) continue;
        // brush square -> circle
        if (dx*dx + dy*dy > brushSize * brushSize) continue;

        int cx = constrain(baseX + dx, 1, N);
        int cy = constrain(baseY + dy, 1, N);

        float amtX = (mouseX - pmouseX) * currentForceAmount;
        float amtY = (mouseY - pmouseY) * currentForceAmount;
        fluid.addVelocity(cx, cy, amtX, amtY);

        float cells = PI * brushSize * brushSize;
        float amountPerCell = currentMilkAmount / cells;
        if (mouseButton == LEFT) {

          fluid.addDensity(cx, cy, amountPerCell);
        } else if (mouseButton == RIGHT) {
          fluid.setDensity(cx, cy, 0);
        }
        if (symmetryEnabled) {
          int mirrorGX = 2 * mugCX - (baseX + dx);   // grid X mirrored
          int mirrorGY = baseY + dy;                 // same Y

          int mx = constrain(mirrorGX, 1, N);
          int my = constrain(mirrorGY, 1, N);
          int mgx = mx - mugCX;
          int mgy = my - mugCY;
          if (mgx * mgx + mgy * mgy <= r2) {

            // invert horizontal component of velocity
            fluid.addVelocity(mx, my, -amtX, amtY);

            if (mouseButton == LEFT) {
              fluid.addDensity(mx, my, amountPerCell);
            } else if (mouseButton == RIGHT) {
              fluid.setDensity(mx, my, 0);
            }
          }
        }
      }
    }
  }
}

public void handleSimulation() {
  drawLayer.beginDraw();
  drawLayer.clear();
  drawLayer.pushMatrix();
  drawLayer.translate(simOffsetX, simOffsetY);

  drawLayer.fill(coffeeColor);
  drawLayer.noStroke();
  drawLayer.rect(0, 0, drawLayer.width, drawLayer.height);

  if (!paused) fluid.step();

  fluid.renderD(drawLayer, SCALE, 98 - renderExposure);

  if (overlayLoaded) {
    drawLayer.imageMode(CORNER);
    drawLayer.image(mugOverlay,
      -simOffsetX,
      -simOffsetY);
  }
  drawLayer.popMatrix();
  drawLayer.endDraw();
}
void draw() {
  handleBrush();
  handleSimulation();

  background(60);
  image(drawLayer, 0, 0);
  fill(255);
  textAlign(LEFT, TOP);
  textSize(12);
  text("FPS: " + nf(frameRate, 0, 1),
    width - 90, height - 30);
}

public void exportPNGButton() {
  File dir = new File(exportFolder);
  if (!dir.exists()) dir.mkdirs();

  String fname = String.format("latte-%04d-%02d-%02d-%02d-%02d-%02d.png",
    year(), month(), day(), hour(), minute(), second());

  drawLayer.save(exportFolder + "/" + fname);
  println("Saved " + fname + "  →  " + exportFolder);
}

public void aboutButton() {
  cp5.getGroup("aboutGroup").show();
}

public void closeAboutBtn() {
  cp5.getGroup("aboutGroup").hide();
}
