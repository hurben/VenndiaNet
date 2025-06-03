Venn-diaNet : venn diagram based network propagation analysis framework for comparing multiple biological experiments
======

Authors
--------

### Benjamin Hur, Dongwon Kang, Sangseon Lee, Ji Hwan Moon, Gung Lee and Sun Kim
https://doi.org/10.1186/s12859-019-3302-7

Key files
--------

### view.py
> Handles I/O between front-end & back-end

### MLV.py 
> Venndianet main program 

1. creates RWR ready format from initial inputs (.mlv)
2. creates adjaceny matrix, p0_matrix after seed is defined by user
3. run RWR.R
4. create visualization format for cytoscape.js

### FL_MLV_Methods_v6.py
> Function Library for MLV

- Currently, Venn-diaNet uses String DB
- If users want to use other network topology, investigate 'line 65: StringDB_to_dict()'

### RWR.R
> Network propagation.


### ```MLV_static```
- all static files are served here. (bootstrap, jquery, fontawesome, img/js/css ...)


### Other sensitive files are not included in the repository
