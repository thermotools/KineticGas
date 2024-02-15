import numpy as np
import matplotlib.pyplot as plt

class ColorGradient:

    def __init__(self, colors):

        if len(colors) < 2:
            raise ValueError('Must supply at least two colors!')

        self.lenc = len(colors[0])
        for c in colors:
            if len(c) != self.lenc:
                raise ValueError('All color tuples must have equal lenght!')

        self.nc = len(colors)
        self.colors = [c for c in colors]
        self.bins = [i / (self.nc - 1) for i in range(self.nc)]
        self.binwidth = 1 / (self.nc - 1)

    def __len__(self):
        return self.nc

    def __call__(self, v):

        if v < 0:
            v = 0
        elif v > 1:
            v = 1

        if v == 0:
            return self.colors[0]
        elif v == 1:
            return self.colors[-1]

        else:
            for i in range(self.nc):
                if v < self.bins[i]:
                    c1 = self.colors[i - 1]
                    lim1 = self.bins[i - 1]
                    c2 = self.colors[i]
                    break

        return tuple(((c2[i] - c1[i]) / self.binwidth) * (v - lim1) + c1[i] for i in range(self.lenc))

class Colormap2D:

    def __init__(self, gradients):
        if len(gradients) < 2:
            raise ValueError('Must supply at least two gradients!')

        self.leng = len(gradients[0])
        for g in gradients:
            if len(g) != self.leng:
                raise ValueError('All color tuples must have equal lenght!')

        self.ng = len(gradients)
        self.gradients = [g for g in gradients]
        self.bins = [i / (self.ng - 1) for i in range(self.ng)]
        self.binwidth = 1 / (self.ng - 1)

    def display(self):
        x_lst = np.linspace(0, 1, 50)
        y_lst = np.linspace(0, 1, 50)
        c_lst = []
        for i, x in enumerate(x_lst):
            c_lst.append([])
            for j, y in enumerate(y_lst):
                c_lst[i].append(self(x, y))
        for i, x in enumerate(x_lst):
            plt.scatter(np.ones_like(y_lst) * x, y_lst, color=c_lst[i], s=40)
        plt.show()

    def __call__(self, v1, v2):

        if v2 < 0:
            v2 = 0
        if v2 > 1:
            v2 = 1

        if v2 == 0:
            g1 = self.gradients[0]
            g2 = self.gradients[0]
            norm_v2 = 0
        elif v2 == 1:
            g1 = self.gradients[-1]
            g2 = self.gradients[-1]
            norm_v2 = 1
        else:
            for i in range(self.ng):
                if v2 < self.bins[i]:
                    g1 = self.gradients[i - 1]
                    g2 = self.gradients[i]
                    norm_v2 = (v2 / self.binwidth) - self.bins[i - 1] / (self.bins[i] - self.bins[i - 1])
                    break


        return ColorGradient([g1(v1), g2(v1)])(norm_v2)