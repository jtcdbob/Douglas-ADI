% CS 5220 Final Project - Douglas-ADI
% Bob Chen
% 2015-11-03

# Introduction

## Douglas Alternating Direction Implicit

The [Alternating Direction Implicit][fw-wiki](ADI) method is a finite difference method for solving parabolic, hyperbolic and elliptic partial differential equations. We use a refined method by Jim Douglas, Jr that is introduced in the paper *Alternating Direction Methods for Three Space Variables*.

[fw-wiki]:https://en.wikipedia.org/wiki/Alternating_direction_implicit_method

## Our Goal

We will start from the current implementation of Douglas API and attempt three tasks:

1.  *Profiling*:  We will find out the bottlenecks in the current code. The goal is to determine what parts of the code are slowest and highlight what we can do to improve the performance in limited amount of time. The existing code has a timer that reflects the total computation time but does not include more details. We will use profiling tools, such as VTune Amplifier (if it becomes available again), and also manually instrument the code with timers.

2.  *Parallelization*: The current code has an attempted parallelization with OpenMP but it has too much overhead to be of any practical use. We will improve the serial code and parallelize it to achieve better performance. Also, we will try to run the parallel code on both the main cores on the nodes and on the Xeon Phi boards for very large problem size. We will set up both strong and weak scaling studies by varying the number of threads/processes we employ.

3.  *Tuning*:  In the end, we will tune the code to get it to run as fast as possible. And we will compare our results to other sparse solvers such as conjugate gradient if time permitted. As the roofline model suggests, we will run into memory issues as we keep increasing the problem size, we will tune the code to deal with this issue, or at least we will propose plans to further scale the code for larger problems.

## Reference
1. Alternating direction implicit method on Wikipedia[fw-wiki].
2. *Alternating Direction Methods for Three Space Variables*[fw-paper].
3. Teaching archive from Prof. Chris Anderson.[fw-archive]

[fw-paper]:http://download.springer.com/static/pdf/464/art%253A10.1007%252FBF01386295.pdf?originUrl=http%3A%2F%2Flink.springer.com%2Farticle%2F10.1007%2FBF01386295&token2=exp=1446587258~acl=%2Fstatic%2Fpdf%2F464%2Fart%25253A10.1007%25252FBF01386295.pdf%3ForiginUrl%3Dhttp%253A%252F%252Flink.springer.com%252Farticle%252F10.1007%252FBF01386295*~hmac=846042c5c990898cd81c5a209b987d766233b4c5ae6ca29b72684e4594ee96be
[fw-archive]:http://www.math.ucla.edu/~anderson/TeachingArchive/index.html
