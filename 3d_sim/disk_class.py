import numpy as np
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from matplotlib.patches import Circle, PathPatch
import mpl_toolkits.mplot3d.art3d as art3d

parser = argparse.ArgumentParser(
        description="""class for 3d protein simulation.
each protein location is defined a vector of its radius, height two angles in respect to center of sphere,
and distance from center of sphere.""",
)

args = parser.parse_args()

class protein:
    def __init__(self, angle1, angle2, dist, Rdisk=7, height=2, center=np.array([0,0,0])):
        # general, same for all disks
        self.Rdisk  = Rdisk     # nm
        self.height = height    # nm
        self.center = center    # 3D coord

        # specific, spherical coordinates 
        self.angle1 = angle1    # rad, theta
        self.angle2 = angle2    # rad, phi
        self.dist   = dist      # nm, distance from center

    def plot_one(self): #Rdis=self.Rdisk, height=self.height, center=self.center):
        # TODO: implement
        #   2. plot disk around localized center, radius R
        #   3. plot two more disks at +0.5h and -0.5h to show disk
        """
        plots one protein disk in 3D space
        """
        ax = plt.figure().add_subplot(projection='3d')
        x = self.dist * np.sin(self.angle1) * np.cos(self.angle2)
        y = self.dist * np.sin(self.angle1) * np.sin(self.angle2)
        z = self.dist * np.cos(self.angle1)

        # plot center of disk
        ax.scatter(xs=0, ys=0, zs=0, color='k')
        ax.scatter(xs=x, ys=y, zs=z, color='b')

        # plot middle disk
        theta = np.linspace(0, 2 * np.pi, num=50)
        x = self.Rdisk * np.cos(theta)
        y = self.Rdisk * np.sin(theta)
        ax.plot(x, y, gapcolor='b', alpha=0.3)
        

        # plot top (+0.5h) and bootom (-0.5h) disks

        plt.show()

def plot_original_sphere(ax, R=30):
    # TODO: implement
    #       source: https://stackoverflow.com/questions/11140163/plotting-a-3d-cube-a-sphere-and-a-vector

    # draw sphere
    u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
    x = R * np.cos(u) * np.sin(v)
    y = R * np.sin(u) * np.sin(v)
    z = R * np.cos(v)
    # alpha controls opacity
    ax.plot_surface(x, y, z, color="g", alpha=0.2)
    return ax

def plot_n_objs():
    # TODO: implement
    return "not implemented yet"

def plot_sim():
    # TODO: implement
    # plots all objects and full sphere
    # ax = plt.figure().add_subplot(projection='3d')
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')


    ax = plot_original_sphere(ax)

    # center of sphere
    ax.scatter(xs=0, ys=0, zs=0, color='k')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_xlabel('Z')
    plt.show()
    return "not implemented yet"

if __name__ == "__main__":
    # TODO: data structure to hold n-objs as a matrix
    obj1 = protein(angle1=np.pi/2, angle2=np.pi/2, dist=20)
    print(f"radius: {obj1.Rdisk} nm")
    print(f"height: {obj1.height} nm")

    obj1.plot_one()
    # plot_sim()
