from __future__ import division
import numpy as np
import os
import sys

#s = 7000
#y_bounds = [,[2.5, 3],[3, 3.5],[3.5, 4],[4, 4.5],[0, 0.75],[0.75, 1.5],[1.5, 2],[2, 2.4],[0, 0.9],[0.9, 1.2],[1.2, 1.6],[1.6, 2.1],[2.1, 2.4]]
nbin = 200
ptmin = 0
ptmax = 20

def desc(s, y_bound, nbin):
    return "# experimental data d^2sigma/dpT/dy\n# energy ylow yup nbin\n{:.1f} {:.1f} {:.1f} {}\n# ptlow ptup value err\n".format(s, y_bound[0], y_bound[1], nbin)

def line(pt):
    return "{} {} 0.0 0.0\n".format(pt, pt+ptmax/(nbin-1))


with open("pseudodata_alice.dat", "w+") as file:
    file.write(desc(5020., [-0.9, 0.9], nbin))
    for pt in np.linspace(ptmin, ptmax, num = nbin):
        file.write(line(pt))
    file.write(desc(5020., [-1.4, 0.4], nbin))
    for pt in np.linspace(ptmin, ptmax, num = nbin):
        file.write(line(pt))
    file.write(desc(8016., [-4.5, -3.0], nbin))
    for pt in np.linspace(ptmin, ptmax, num = nbin):
        file.write(line(pt))
    file.write(desc(8016., [2.0, 3.5], nbin))
    for pt in np.linspace(ptmin, ptmax, num = nbin):
        file.write(line(pt))






"""
# experimental data d^2sigma/dpT/dy
# energy ylow yup nbin
7000. 0. 0.9 10
# ptlow ptup value err
8. 9. 42.83305227655986 2.2936712493651674
9. 10. 26.30691399662732 1.8929126745905267
10. 11. 16.188870151770658 1.011804384485666
11. 12. 10.320404721753793 0.6184526415140724
12. 13.5 5.919055649241146 0.33768945016021557
13.5 15. 3.102866779089376 0.159980320075972
15. 18. 1.4182124789207418 0.0746003249935739
18. 30. 0.23440134907251262 0.013595712897636678
30. 45. 0.013153456998313658 0.0008762483006250645
45. 70. 0.0011804384485666105 0.00031413045548393185
"""