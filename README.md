# bezierplot
Given a smooth function, bezierplot returns a smooth bezier path written in [tikz](https://www.ctan.org/pkg/pgf). It finds special
points such as extreme points and inflection points and reduces the number of used points.

![bezierplot-nodes1](https://user-images.githubusercontent.com/11213578/37959146-bc034794-31b2-11e8-8f90-5726c42cf6c9.png)

The upper graph of sqrt(x) used bezierplot, the lower used the built-in
plotting function of tikz with 100 samples (no smoothing) and is still
quite inexact at the beginning.

![bezierplot-nodes2](https://user-images.githubusercontent.com/11213578/37959147-bc2429fa-31b2-11e8-82ae-a988700cd7f0.png)

For more information read the documentation.

## Author

Linus Romer
