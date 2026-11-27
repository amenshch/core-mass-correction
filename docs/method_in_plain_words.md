# The correction method in plain words

A guide to what the method does, what the borrowed vocabulary means, and where
the new column-density calibration would fit. Nothing here is new physics; it is
the same procedure described without jargon.

---

## 1. The problem, in one paragraph

We extract a core from a Herschel map, fit a modified blackbody to its fluxes,
and obtain a mass. That mass is wrong, for reasons we understand: part of the
core lies outside the measurement footprint, the background subtracted under the
core is not the right background, and a single temperature cannot describe an
object with a cold interior and a warm skin. We cannot measure the error,
because for a real core nobody knows the true mass. So we build models where the
true mass *is* known, put them through the same measuring procedure, see how
wrong the answer comes out, and use that to correct real cores.

---

## 2. What "training" actually means

The word sounds like something the computer does on its own. It is much simpler
than that.

**Training is building a table.** Each row of the table is one model, measured in
the same way a real core is measured. A row contains:

- the few numbers we could also measure on a real source — call these the
  *observables*;
- the answer we want but cannot measure on a real source — for the mass
  correction, the factor by which the measured mass must be multiplied.

That is the whole of it. Our table has about 30 000 rows, built from 1597
models, each model measured through many different footprint sizes.

**Using the table** is equally simple. A real core arrives with its own
observables. We look through the table for the rows whose observables are closest
to the real core's, and read off their answers. If the closest models needed to
be multiplied by 1.3, the real core probably needs about 1.3 too.

There is no learning, no intelligence, no black box. It is a lookup table with
interpolation. Everything else is detail about how the lookup is done.

**Two refinements, and why they matter.**

*Interpolation instead of nearest-neighbour.* Rather than taking the single
closest row, we take the 320 closest and fit a straight-line trend through them,
then evaluate that trend at the real core's position. This is called *local
linear regression*: "local" because only nearby rows are used, "linear" because a
straight line (in four dimensions, a flat plane) is fitted through them, and
"regression" is simply the statistician's word for fitting. The benefit over
picking the nearest row is that the answer varies smoothly, and the fit can
follow a trend between grid points instead of stepping from one to the next.

*Distances must be fair.* The observables have different units and different
spreads, so before measuring "closeness" each one is divided by its own spread
across the table. Otherwise the observable with the largest numbers would decide
which rows count as near.

---

## 3. The observables, and why these four

For the lookup to work, the observables must satisfy three conditions: they must
be measurable identically on a model and on a real source; they must not depend
on distance, since the models sit at one distance and the clouds at seven; and
together they must actually determine the answer.

The four in use are:

| Observable | Meaning |
|---|---|
| background column density | how bright the cloud is under the source |
| concentration | the source's peak divided by its mean inside the footprint |
| footprint factor | the footprint size divided by the source's own width |
| footprint radius | the angular size of the footprint, in arcseconds |

The middle two are ratios, so they carry no units and no distance. The last is an
angle, not a physical length, and that is deliberate: using a physical size in
parsecs introduced a false trend with distance across the seven clouds, because
the same physical size looks different at 130 and 1700 pc.

The fourth observable is the one that does most of the work. Without it the
correction is systematically off by about 13 per cent and cannot handle the
models whose outer edge falls below the surrounding cloud; with it, the
correction is unbiased and every model is recovered within a factor of two.

---

## 4. Why we predict a *ratio*, never a quantity

This is the single most useful idea in the method, and it generalises.

If we ask the table for the *mass* of a real core, it does badly: an error of
about 30 per cent. If we ask it for the *factor* by which the measured mass must
be multiplied, it does very well: about 6 per cent.

The reason is that the measured mass already contains most of the information.
The measured flux, the measured size, the fitted temperature — all of that is
real information about this particular source. Predicting the mass from scratch
throws it away and relies on the models alone. Predicting a correction factor
keeps it and asks the models only for what the measurement got wrong, which is a
much smaller and much smoother thing.

The same held for every other quantity we tried. The stability parameter is
recovered twice as well as a ratio than directly. The mass-averaged temperature,
twice as well. Where a ratio was not possible, accuracy was poor.

**Rule of thumb:** ask the models for a correction, not for an answer.

---

## 5. How we know it works: leave-one-out

An obvious objection: of course the table reproduces its own rows. So the test
is done with the answer hidden.

Take one model out of the table entirely — all of its rows. Rebuild the lookup
from the remaining 1596 models. Now ask it about the model we removed, and
compare with that model's known truth. Repeat for every model in turn.

This is *leave-one-out*. It measures what happens to a source the table has never
seen, which is the situation of every real core. Our result: the corrected mass
is unbiased, with a typical error of 6 per cent and every model within a factor
of two.

A second, independent check used a separate grid of 1400 models built at
parameter values deliberately placed *between* those of the training models.
Those were never in the table at all, and gave the same accuracy. So the
interpolation between grid points is sound.

**An important honesty point.** On real extractions, the accuracy is not 6 per
cent but about 30 per cent. The extra error is not in the correction; it is in
the extraction — noise, blending, imperfect background subtraction. The method is
five times more precise than the data it is applied to. This means further
refinement of the correction is pointless; the limit lies elsewhere.

---

## 6. The chain as it runs today

```
Herschel maps (160-500 um)
        |
        v
  hires  -->  column-density map at 13.5"
        |
        v
  getsf  -->  source catalogue: positions, fluxes, sizes, footprints,
        |     background under each source
        v
  the pipeline reads the catalogue and, for each source:
        |
        +--> fits a modified blackbody to the fluxes  -->  measured mass
        |
        +--> forms the four observables from the catalogue
        |
        +--> looks them up in the model table          -->  correction factor
        |
        v
  corrected mass = measured mass x correction factor
```

Everything from "the pipeline reads the catalogue" downward is in
`mass_correction_pipeline_6.py`. The model table is built once, at start-up, from
two files: the model catalogue and the recovery tables.

---

## 7. Where the new column-density correction fits

**It is external, and it is upstream of everything.** It is not in the pipeline
and cannot be, because it acts on the *map*, before getsf ever runs.

The finding behind it: the column-density map that hires produces is
systematically too low, because the same single-temperature assumption that
biases core masses also biases every pixel. The size of the bias depends on how
much material is along the line of sight — about 8 per cent at 7.5x10^20 cm^-2,
21 per cent at 6x10^21, and 52 per cent at 4.8x10^22 — and it is worse again at
the centre of a core that the beam does not resolve.

**The correction works the same way as the mass correction**, which is why it is
worth understanding the vocabulary once. A table is built from model stamps. Each
row is one pixel. The observables are quantities hires itself produces. The
answer is the factor by which that pixel's column density must be multiplied.

The observables are the interesting part. Hires internally fits the dust
emission several times over, at different resolutions and with different band
combinations, and then discards those intermediate results. They are exactly what
is needed, because **where the fits disagree with each other, the line of sight is
not isothermal or not resolved — and that is precisely where the bias is large.**
Three are used:

- the final hires column density;
- the ratio of the column fitted at 18.2" (160-250 um) to that fitted at 36.3"
  (all four bands) — this measures beam dilution of an unresolved peak;
- the difference of the corresponding temperatures.

Accuracy, tested by leaving one model out at a time:

| where | error before | error after |
|---|---|---|
| cloud, away from sources | 21 % | 0.4 % |
| under a typical Aquila core | 21 % | 0.8 % |
| brighter cores | 27 % | 1.6 % |
| very bright, compact peaks | 50 % + | 6-9 % |

---

## 8. What would have to be produced for it to work

**On the hires side, one change:** keep the intermediate fits instead of deleting
them. For each map we would need, alongside the final column-density image:

- column and temperature fitted at 18.2" from 160 and 250 um;
- column and temperature fitted at 36.3" from all four bands;
- optionally the same at 36.3" from 350 and 500 um, which adds a smaller further
  improvement.

**On the analysis side, two new files:**

1. *The calibration table* — the equivalent of the model catalogue and recovery
   tables, but for pixels: one row per pixel of the model stamps, holding the
   three observables and the required correction factor. Built once. From 16
   models it has about 650 rows; more models would make it more robust,
   especially for compact bright peaks, where there are currently only 17 rows.

2. *A script that applies it* — reads the calibration table, reads the four hires
   images of a real field, evaluates the correction at every pixel, and writes a
   corrected column-density image. The lookup code already exists and is the same
   one the mass correction uses.

**The order of operations would become:**

```
Herschel maps --> hires --> corrected column-density map --> getsf --> pipeline
```

with everything downstream unchanged.

**One further consequence, and it is a benefit.** Today the model table is built
from the models' *true* column densities, while real sources are measured on a
hires map that is biased. The two are not the same quantity, and we measured what
that costs: a few per cent in the median correction factor, up to 13 per cent in
the worst cases. If real maps were corrected first, the two would be the same
kind of quantity by construction, and that mismatch would disappear.

---

## 9. Vocabulary, once and for all

| Term | Plain meaning |
|---|---|
| training | building the table of models with known answers |
| training set | the table itself |
| observables / features / axes | the few numbers used to find similar models |
| target | the answer being looked up |
| local linear regression | fit a straight-line trend through the nearest 320 rows |
| neighbours | the rows closest in observable space |
| leave-one-out | hide one model, rebuild the table, test on it, repeat |
| out of training | the source being corrected is not itself in the table |
| interpolation | reading between the grid points |
| extrapolation | asking outside the range the table covers — unreliable |
| convex hull / domain | the region of observable space the table covers |
| bias | the systematic part of an error, the median offset |
| scatter | the random part, how much individual answers vary |
| transfer | whether something that works on models also works on real data |

---

## 10. Five things worth remembering

1. Training is building a table of models; using it is looking up similar models
   and interpolating. There is nothing else in the machinery.

2. Ask the models for a correction factor, never for a quantity. The measurement
   already contains most of the information.

3. Test by hiding the answer. Anything else measures how well the table
   reproduces itself.

4. What works on models need not work on real data. Twice in this project a large
   model-side gain vanished on real extractions, because the same quantity was
   measured slightly differently on the two sides. Small differences in *how* a
   number is produced matter more than large differences in the number itself.

5. The real-data accuracy of the mass correction is set by the extraction, not by
   the correction. The place to improve is upstream — which is exactly what the
   column-density calibration would be.
