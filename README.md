# BIOL3207 Species Distribution Modelling Workshop (2026)

Workshop materials for the two-week Species Distribution Modelling (SDM)
module of **BIOL3207 – Evolution, Ecology and Genetics** at the Australian
National University.

The workshops walk students through a full SDM pipeline using the frilled
lizard (*Chlamydosaurus kingii*) as a case study:

- **Week 7 — Species occurrences & geometric SDMs.** Downloading and
  cleaning GBIF occurrence data, visualising points and polygons, and
  fitting simple geometric range models (minimum convex polygons, alpha
  hulls).
- **Week 8 — Correlative SDMs & projections.** Generating background /
  pseudo-absence points, fitting correlative models with environmental
  predictors, validating them, and projecting habitat suitability through
  time under climate change scenarios.

## Repository contents

| File / folder | What it is |
| --- | --- |
| `workshop_instructions_wk7.qmd` / `.html` | Week 7 instructions (Quarto source + self-contained HTML) |
| `workshop_instructions_wk8.qmd` / `.html` | Week 8 instructions (Quarto source + self-contained HTML) |
| `data/` | Raw input data for both workshops (occurrences, climate layers, shapefiles) |
| `style.css` | Custom styling applied to the rendered HTML |
| `frilled_lizard.jpg` | Cover image used in the instructions |

## Getting started (students)

1.  **Clone or download this repository.**
    - *GitHub Desktop:* File → Clone Repository → URL →
      `https://github.com/alexskeels/BIOL3207_SDM_Workshop`
    - *Or* click the green **Code → Download ZIP** button on GitHub and
      unzip locally.
2.  **Create an R Project** in the cloned folder (*File → New Project →
    Existing Directory* in RStudio). Opening the `.Rproj` automatically
    sets your working directory to the repo root, so every relative path
    used in the workshop (e.g. `read.csv("data/occurrences.csv")`) will
    just work.
3.  **Check the `data/` folder is populated.** If it's empty or missing
    files, download `data.zip` from the corresponding Canvas module and
    unzip its contents into `data/`.
4.  Open `workshop_instructions_wk7.html` (or `wk8`) in a browser and
    follow along. Full package install instructions are in the
    *Packages for today* section of each document.

## Rendering the instructions

The `.qmd` files use `embed-resources: true`, so a single `quarto render`
produces a standalone `.html` with all CSS, JS, fonts, and images baked
in — suitable for uploading to Canvas as a File.

```bash
quarto render workshop_instructions_wk7.qmd
quarto render workshop_instructions_wk8.qmd
```

## Licence & attribution

Materials © Alex Skeels, for use in BIOL3207 at ANU. Please get in touch
before reusing in other teaching contexts.
