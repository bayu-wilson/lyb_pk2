#!/usr/bin/env python

import imageio

nqso = 100

images = []
for qidx in range(nqso):
    filename = "sightlines/q{:02}.png".format(qidx)
    images.append(imageio.imread(filename))
imageio.mimsave('figures/all.gif',images,fps=15)



#
# import imageio
# im1 = "IMAGE_1.jpeg"
# x = imageio.imread(im1)
# images = []
# for i in range(275):
# 	ADD= 10 * i
# 	images.append(x[0+ADD:250+ADD,0+ADD:250+ADD])
# imageio.mimsave('movie.gif', images, duration = 0.1)

"""
filenames = [im1,im2]
images = []
for filename in filenames:
    images.append(imageio.imread(filename))
imageio.mimsave('movie.gif', images)


import imageio
images = []
for filename in filenames:
    images.append(imageio.imread(filename))
imageio.mimsave('/path/to/movie.gif', images)
"""
