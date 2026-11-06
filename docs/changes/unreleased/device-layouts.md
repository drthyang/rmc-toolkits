### Layouts for iPhone and iPad

- Phones and portrait iPads: the run-folder field and the bond-window and bin-width fields now
  use 16px text, so iOS no longer zooms the page when one takes focus (a rule meant to do this
  lost to the fields' own sizes).
- The Copilot's sheet keeps its composer on screen on a phone in landscape. An open sheet takes
  at least 19rem (up to 72% of the window); before, 46% of a ~400px-tall window cut off the
  composer under the sheet's head.
- Displacement Directions on a touch screen: a tap on the sphere reads out that cell (direction,
  enhancement, atoms, local z), as the mouse pointer does. A tap off the sphere or a drag clears
  it, and past the middle of the view the readout sits left of the point, so it stays on the
  card.
- iPad 11" in portrait: PCA Ellipsoid shows its two 3D views side by side, under statistics
  laid out 2 × 2 (also on a 13" in portrait), where one long column used to put the views a
  screen down. In landscape, the four statistics columns are sized so the crystal orientation
  no longer scrolls sideways.
- iPad 11" in landscape and 13" in portrait: Displacement Directions shows the sphere and the site
  picker side by side, with the axis views in a row under both.
- Windows 1041–1319px wide (an iPad 11" in landscape, beside the Copilot rail): the header takes
  two tidy rows, brand and data source then the tabs. Before, the data source wrapped alone onto
  the right of a second row.
- Touch screens: tap targets are 44px. Fields, selects and header buttons are 44px tall, and the
  page tabs too. Small buttons keep their look but respond to taps across a 44px strip: tool and
  save buttons, the 3D and card toggles, pills, the Copilot's buttons, and the footer links.
  Sliders respond across a 44px-tall strip, and the space-group ladder's bricks are 44px.
