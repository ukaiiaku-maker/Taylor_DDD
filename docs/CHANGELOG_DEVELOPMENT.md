# Development changelog summary

## v12

Pinned lines retain raw curvature/backstress mobility by default. This addressed the issue where center-of-mass projection could remove the contact-loading force and suppress Taylor hardening.

## v13

Plastic strain was corrected to use actual swept line motion after capture, depinning, snapping, and relaxation. This fixed a bookkeeping problem where FAC changes altered event counts without correctly updating stress/strain feedback.

## v14

Diagnostics were expanded: before-step and after-step stress, step stress increments, plastic ratio, preflight barrier/rate diagnostics, and run summaries. A field-name and after-step stress write-order issue was fixed in the checked v14 version.

## v15

Branch-specific barrier floors were added for forest crossing and Peierls branches. This improved event resolution but revealed that many high-density depin events were still local-stress-cap dominated.

## v16

Avalanche/burst diagnostics were introduced: burst grouping, CCDF tables, randomized null comparison, and tail-fit comparisons.

## v17

The active local stress cap was removed and the local barrier stress was computed using a physical feeding length. This made the model cap-free in the relevant diagnostics and made the high-density hardening/intermittency analysis physically interpretable.
