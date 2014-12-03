================================================
SOURCE README FOR SENSITIVE COUTURE
Nobuyuki Umetani (umetani@ui.is.s.u-tokyo.ac.jp)
================================================


REFERENCES
==========
This code is the implementation of the following paper. Note that implementation may not be exact. Please cite this paper if you use our implementation.

N.Umetani, D.Kaufman, T.Igarashi, E.Grinspun, "Sensitive Couture for Interactive Garment Editing and Modeling", Transaction of ACM SIGGRAPH, Vol. 30, No. 4, 2011



Required 
========
We have been running our demo using Core 2 Duo 2.66 GHz 64bit machine with 4Gb RAM. We use OpenGL in the program, so graphical hardware is necessary.


Compiling
=========
The command "make" should compile and link the program.


The Design Window
===================

Tools
-----
To change the tool, press right button and select tool from pop-up menu.
+ Drag : you can drag vertex/edge/loop of a pattern.
+ Dart Cutter : you can draw dart on a pattern. 
+ Edit Curve : you can grab the curve by click and dragging a polyline edge.
+ ChangeToLine : you can change a polyline edge to line.
+ ChangeToPolyline : you can change a line edge to a polyline edge.
+ SmoothPolyline : you can smooth a polyline edge by moving mouse with left button pressed near the polyline edge.

Keyboard
------
+ ' '  : change 3D model
+ 'b'  : Initialize the simulation 

View Change
-----------
This design window should be activated so as to change the view of the design window.
+ Pan : Shift+Left drag
+ Zoom : Page Up / Down



The Simulation Window
=====================

Texture
-------
you can change texture of cloth from pop up menu. select [texture] and chose one from submenu.

View Change
------
This design window should be activated so as to change the view of the simulation window.
+ Pan : Shift+Left drag
+ Rotataion : Ctrl+Left drag



