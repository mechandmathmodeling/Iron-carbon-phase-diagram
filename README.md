# Iron-carbon phase diagrams

## The one-sentence version

Iron can hold a little carbon in its crystal gaps, and how much it can hold
depends on temperature. The diagram is a map of that: at a given temperature
and carbon content, it tells you what the material is made of.

## Why carbon matters at all

Imagine iron atoms as a stack of oranges, and carbon as a marble small enough
to slip into the gaps between oranges. Pure iron has two stacking patterns.
At low temperature the oranges pack in a loose cubic arrangement (BCC, called
ferrite, α) with small gaps. At higher temperature they pack in a tighter
arrangement (FCC, called austenite, γ) with bigger gaps. The bigger gaps hold
more carbon marbles.

That is the whole story of steel. Austenite holds up to 2.14 wt% C. Ferrite
holds only 0.022 wt% C at its most. Whatever carbon doesn't fit has to go
somewhere else, and "somewhere else" is a separate compound, cementite (Fe₃C),
which is 6.67 wt% carbon.

## Reading the diagram

Put temperature on the vertical axis and carbon content on the horizontal
axis. Each line is a boundary: cross it and the material changes what it is
made of. The lines are not decoration; they are the answer to the question
"how much carbon can each phase hold at this temperature?"

Ask this question of any point you pick. The region you land in tells you the
phases. The lines between regions tell you the composition of each phase.

## The lever rule is a seesaw

If you land in a region with two phases, you need to know how much of each
there is. Imagine a horizontal line (a tie line) drawn across the two-phase
region at your temperature. Its ends sit on the boundaries, and they tell you
the composition of each phase. Your alloy's composition, C₀, sits somewhere
on that line.

Now think of the line as a seesaw balanced on the fulcrum at C₀. Each phase
is a weight, and its weight is its mass fraction. The rule is simply that the
seesaw balances:

    fraction_α × (C₀ − C_α) = fraction_γ × (C_γ − C₀)

The fraction of a phase is decided by the arm on the *opposite* side. If your
alloy is close to C_α, the arm to C_α is short, so the α phase must be heavy
to balance. That is why a carbon-poor alloy is mostly the phase that holds
little carbon.

Worked example: alloy at 0.40 wt% C just below 727 °C. The tie line runs from
ferrite (0.022) to austenite (0.76). Left arm 0.378, right arm 0.36 in
composition units, so the fractions are 0.378/0.738 ≈ 51% austenite
(which becomes pearlite) and 49% ferrite. The app shows this as a balanced
seesaw with the torque numbers printed under it.

## The dot splits in two

As the alloy cools, its point moves straight down the vertical line at C₀.
While it stays in a single-phase region, it is one dot. When it enters a
two-phase region, the dot is no longer one thing. It becomes two dots, one on
each boundary, linked by the tie line. The two dots tell you the composition
of each phase, and the seesaw tells you how much of each there is.

That is the picture to hold in your head: one alloy, two phases, two
compositions, one balance.

## What happens below 727 °C

This is the part people get wrong, so think it through slowly.

The eutectoid reaction at 727 °C is where austenite (0.76 wt% C) can no
longer hold its carbon. It splits into ferrite and cementite. Below 727 °C
there is no austenite left. The only phases are ferrite and cementite.

Here is the key point. This is an equilibrium diagram, so the compositions
keep changing as the temperature drops:

- **Cementite** is a compound with a fixed composition, 6.67 wt% C. It never
  changes.
- **Ferrite** is a solid solution, and its carbon solubility follows the
  solvus curve. It is 0.022 wt% C at 727 °C and drops to almost zero at room
  temperature. As it cools, ferrite rejects carbon, and that carbon goes into
  cementite.

So the lever rule below 727 °C uses a moving endpoint, C_α(T):

    fraction_Fe3C = (C₀ − C_α(T)) / (6.67 − C_α(T))

For 0.40 wt% C, that gives about 5.7% cementite at 727 °C and about 5.9% at
room temperature. The change is small, but it is real, and it comes from the
solvus.

## Pearlite is a recipe, not a phase

Pearlite is a layered mixture of ferrite and cementite. It forms when
austenite of eutectoid composition transforms at 727 °C. It is a
microconstituent, a picture of how the structure looks under a microscope.
It is not a phase on the diagram, and it has no composition of its own.

To describe the microstructure of a hypoeutectoid steel (C < 0.76 wt%)
after cooling past 727 °C, you split the material into two parts: the
pearlite made from the austenite, and the proeutectoid ferrite that formed
before the eutectoid reaction.

    fraction_pearlite = (C₀ − 0.022) / (0.76 − 0.022)

Hypereutectoid steels (0.76 < C < 2.14 wt%) work the same way, with
proeutectoid cementite:

    fraction_proeutectoid_Fe3C = (C₀ − 0.76) / (6.67 − 0.76)

The app draws the first form of the picture from these numbers. It is the
727 °C reference picture. In slow equilibrium cooling, the ferrite inside
pearlite gives up a little carbon, which the app does not model.

## Things to check before you trust the app

- **The boundary curves are approximations.** The liquidus, solidus, GS (A3),
  SE (Acm), and solvus arrays are hand-placed. Use a reference diagram for any
  number that gets graded.
- **δ-ferrite and the peritectic are not modeled.** Above about 1394 °C, low-
  carbon alloys actually contain δ-ferrite, which the app treats as austenite.
- **Liquid-phase details are simplified.** The liquid + austenite region uses
  rough liquidus and solidus curves.
- **The microstructure picture is a schematic.** Grains are random cells
  colored by phase fraction. Pearlite is stripes, not true lamellar spacing.
- **Steels only.** The composition slider stops at 2.14 wt% C.

## Check yourself

1. Why is austenite able to dissolve more carbon than ferrite?
2. An alloy at 0.20 wt% C sits just below 727 °C. Which phases are present,
   and what is the pearlite fraction in its microstructure?
3. Why does the cementite fraction below 727 °C change with temperature, even
   though cementite's own composition never changes?
4. Move C₀ along a tie line. Which arm gets longer, and which phase gets
   heavier?

If you can answer these in plain words, without looking, you understand the
diagram. The seesaw is the lever rule, the dot splitting is the tie line, and
pearlite is the recipe that describes how the phases were arranged, not a
phase of its own.