
import jax.numpy as jnp
from jax import grad, jit
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

curveFunction = lambda t: jnp.vstack((jnp.cos(
    2 * jnp.pi * t), jnp.sin(2 * jnp.pi * t), 0.3 * jnp.sin(4 * jnp.pi * t)))
# p1 = [2, 0, 0]
# p2 =[18, 0, 0]
# a = 5
# curveFunction = lambda t: jnp.vstack((
#     p1[0] + t * (p2[0] - p1[0]),  # x-coordinate
#     jnp.zeros_like(t),            # y-coordinate, assumed to be 0
#     p1[2] + a * jnp.cosh((p1[0] + t * (p2[0] - p1[0]) - (p1[0] + p2[0]) / 2) / a)  # z-coordinate
# ))

class Curve:

    def __init__(self,
                 curveFunction,
                 nSample=100,
                 bendModulus=1,
                 twistModulus=1,
                 totalTwist=50):
        self.bendModulus = bendModulus
        self.twistModulus = twistModulus
        self.totalTwist = totalTwist
        self.nSamples = nSample
        self.stiffness = 100
        self.dt = 1e-4
        t = jnp.linspace(0, 1, self.nSamples + 1)
        self.verts = curveFunction(t)[:, :self.nSamples].reshape(
            3, self.nSamples)

        self.velocities = jnp.zeros((3, nSample))
        self.rst = self.calculate_initial_lengths()

    def calculate_initial_lengths(self):
        edges = jnp.roll(self.verts, -1, axis=1) - self.verts
        lengths = jnp.sqrt(jnp.sum(edges**2, axis=0))
        return lengths

    def computeEnergy(self, verts):

        edges = jnp.roll(verts, -1, axis=1) - verts
        lengths = jnp.sqrt(jnp.sum(edges**2, axis=0))
        stretch = lengths - self.rst
        elasticEnergy = 0.5 * self.stiffness * jnp.sum(stretch**2)

        edgesR = jnp.roll(verts, -1, axis=1) - verts
        edgeLengthsR = jnp.sqrt(jnp.sum(edgesR**2, axis=0))
        dualLengths = 0.5 * (edgeLengthsR + jnp.roll(edgeLengthsR, 1))
        totalLength = jnp.sum(dualLengths)
        curvatureBinormals = (
            2 * jnp.cross(jnp.roll(edgesR, 1, axis=1), edgesR, axis=0) /
            (jnp.roll(edgeLengthsR, 1) * edgeLengthsR +
             jnp.sum(jnp.roll(edgesR, 1, axis=1) * edgesR, axis=0)))
        bendingEnergy = jnp.sum(self.bendModulus * curvatureBinormals**2 /
                                (2 * dualLengths))
        twistEnergy = self.twistModulus * self.totalTwist**2 / (2 *
                                                                totalLength)
        kineticEnergy = 0.5 * jnp.sum(dualLengths * jnp.sum(self.velocities**2, axis=0))
        return bendingEnergy + twistEnergy + elasticEnergy 


c = Curve(curveFunction)
energy_fn = jit(lambda verts: c.computeEnergy(verts))
totalEnergy_grad = grad(energy_fn, argnums=0)


def update():
    for _ in range(2):
        acceleration = -totalEnergy_grad(c.verts)
        c.velocities += c.dt * acceleration
        c.verts += c.dt * c.velocities

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(c.verts[0, :], c.verts[1, :], c.verts[2, :])

def updates(frame):
    update()
    scatter._offsets3d = (c.verts[0, :], c.verts[1, :], c.verts[2, :])
    return scatter

ani = FuncAnimation(fig, updates, frames=100, blit=False, interval=10)
plt.show()