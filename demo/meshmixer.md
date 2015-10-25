# Extrude a Surface with Meshmixer

How to use Meshmixer's `Select` and `Extrude` tool to solidify a surface mesh created with `phstl`:

1. ![Import Surface](Screenshot-Autodesk Meshmixer - surface.stl.png)
   Open Meshmixer and import [surface.stl](surface.stl)

2. ![Select Surface](meshmixer-select.png)
   Select the `Select` tool and make an interactive selection on the surface. Then choose `Select All` from the `Modify` selection submenu.
   
   ![All selected...](Screenshot-Autodesk Meshmixer - surface.stl-1.png)

3. ![Extrude selection](meshmixer-extrude.png)
   Next choose `Extrude` from the `Edit` selection submenu.

4. ![Extrude settings](Screenshot-Autodesk Meshmixer - surface.stl-2.png)
   Set the offset to a value like `15` mm and set the extrusion `Direction` to `Y Axis` to ensure uniform extrusion.

5. Click `Accept` and export the extruded solid.

# Extrude a Surface with a Flat Base

![Flat base settings](Screenshot-Autodesk Meshmixer - surface.stl-3.png)

Follow the same steps as above, but set the extrusion `EndType` to `Flat`. You may want to specify a negative `Offset` (eg, `-15`) to position the flat base correctly.
