# pyvenn: Venn diagrams for 2, 3, 4, 5, 6 sets

Please refer to the Jupyter notebook for demos and a brief explanation of the
interface; a more complete documentation is in the works as the project keeps
evolving:  
https://github.com/LankyCyril/pyvenn/blob/master/pyvenn-demo.ipynb

This library is an evolution of tctianchi's pyvenn package (see fork URL).  
Their liberal license (Unlicense) allowed me to fork the repository,
change the license to GPLv3, modify the package's interface and, hopefully,
significantly contribute to and improve the library, and make it installable
from PyPI.

The main methods in this version are different from the ones in tctianchi's
implementation, but the original methods are still provided for backwards
compatibility, and I would like to emphasize the importance of tctianchi's work
that allowed for this library to exist (among other things, figuring out the
coordinates best fit for plotting the diagrams' shapes and petals' labels).

This iteration of the library implements two main functions:
* `venn(dataset_dict, **kwargs)` which plots true Venn diagrams for any number
of sets between 2 and 5 using ellipses, and for 6 sets using triangles
* `pseudovenn(dataset_dict, **kwargs)` which plots a Venn-like intersection of
six circles (not all intersections are present in such a plot, but many are).
