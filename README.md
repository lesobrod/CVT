# Counting centroidal Voronoi tesselations

A centroidal Voronoi tessellation ([CVT](https://people.sc.fsu.edu/~jburkardt/classes/urop_2016/du_faber_gunzburger.pdf)) is a special type of Voronoi tessellation in which   
the generating point of each Voronoi cell is also its centroid (center of mass).  
We should highlight _stable_ CVTs (SCVTs) which are local minimizers of energy-like function (see details [here](https://www.microsoft.com/en-us/research/wp-content/uploads/2016/12/On-Centroidal-Voronoi-Tessellation-Energy-Smoothness-and-Fast-Computation.pdf)).  
Unfortunately the [article in wiki](https://en.wikipedia.org/wiki/Centroidal_Voronoi_tessellation) is a stub and contains an extra diagram that is not stable CVT for the square.  

Here I will collect materials related to CVT, and primarily to **counting SCVTs for square and disk**.   
### What are ready now  
- [Tool to detect distinct patterns](DetectPatterns.nb) using Wolfram Mathematica  
(with detailed explanations)
- [Tool for counting](CVTCountingSimple.nb) distinct stable CVTs using Lloyd algorithm  
(as code with short comments)
- [Results](Results) for number of seeds up to 22,  
which are proposed as drafts on OEIS ([disk](https://oeis.org/draft/A363822)  [square](https://oeis.org/draft/A366544))
### TODO list
- Possible problems and pitfalls of generation and counting patterns
- Interactive demonstration of energy function
- Alternatives for initialization and iteration algorithms
- Possible approaches to *unstable* CVT counting
- Realization on Python and Julia
