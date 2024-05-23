# MIT License
#
# Copyright (c) 2024 TriLite https:#github.com/MeshLite/TriLite
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import unittest
import numpy as np
import trilite as TL
import sys


class TestTrimesh(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mesh = TL.Trimesh()

        v0 = np.array([0.0, 0.0, 0.0])
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([1.0, 1.0, 0.0])
        v3 = np.array([0.0, 1.0, 0.0])

        cls.mesh.AddFace([v0, v1, v2])
        cls.mesh.AddFace([2, v3, 0])

    def test_trimesh_functions(self):
        mesh = self.__class__.mesh
        # Basic properties
        self.assertEqual(mesh.NumHalfedges(), 6)
        self.assertEqual(mesh.NumVertices(), 4)
        self.assertEqual(mesh.NumFaces(), 2)

        # Halfedges, Vertices, Faces
        self.assertEqual(len(mesh.Halfedges()), 6)
        self.assertEqual(len(mesh.Vertices()), 4)
        self.assertEqual(len(mesh.Faces()), 2)

        # Positions
        self.assertEqual(len(mesh.Positions()), 4)

        # Halfedge navigation
        self.assertEqual(mesh.HNext(0), 1)
        self.assertEqual(mesh.HPrev(1), 0)
        self.assertEqual(mesh.HStart(3), 2)
        self.assertEqual(mesh.HEnd(3), 3)
        self.assertEqual(mesh.HFace(4), 1)
        self.assertEqual(mesh.HOpposite(5), 2)
        self.assertEqual(mesh.HNextAroundStart(0), 5)
        self.assertEqual(mesh.HPrevAroundStart(5), 0)
        self.assertEqual(mesh.HNextAroundEnd(2), 4)
        self.assertEqual(mesh.HPrevAroundEnd(4), 2)

        # Geometric properties
        self.assertTrue(
            np.array_equal(mesh.HGeometry(0), np.array([1.0, 0.0, 0.0]))
        )

        self.assertEqual(len(mesh.HConnectionsAroundStart(0)), 2)

        self.assertEqual(len(mesh.HHalfedgesAroundHole(3)), 4)
        self.assertAlmostEqual(mesh.HLength(0), 1.0)

        # Vertex properties

        self.assertEqual(mesh.VStarting(0), 5)
        self.assertEqual(mesh.VEnding(0), 4)
        self.assertEqual(len(mesh.VStartings(0)), 2)
        self.assertEqual(len(mesh.VEndings(0)), 2)
        self.assertEqual(len(mesh.VFaces(0)), 2)
        self.assertTrue(
            np.array_equal(mesh.VPosition(0), np.array([0.0, 0.0, 0.0]))
        )
        self.assertTrue(np.allclose(mesh.VNormal(0), [0, 0, 1]))
        self.assertEqual(mesh.VValence(0), 2)
        self.assertTrue(mesh.VIsManifold(0))

        # Face properties
        self.assertEqual(mesh.FHalfedge(0), 0)
        self.assertEqual(len(mesh.FHalfedges(0)), 3)
        self.assertEqual(len(mesh.FVertices(0)), 3)
        self.assertEqual(len(mesh.FPositions(0)), 3)
        self.assertTrue(np.allclose(mesh.FNormal(0), [0, 0, 1]))
        self.assertAlmostEqual(mesh.FArea(0), 0.5)

        # Edge properties
        self.assertEqual(len(mesh.EdgeHalfedges(0)), 1)
        self.assertEqual(len(mesh.EdgeFaces(0)), 1)
        self.assertTrue(mesh.EdgeIsManifold(0))

        # Boundary halfedges
        self.assertGreaterEqual(len(mesh.BoundaryHalfedges()), 0)

        # Modify mesh
        mesh.AddFace(
            [
                np.array([0.0, 1.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([1.0, 1.0, 0.0]),
            ]
        )

        self.assertEqual(len(mesh.RemoveFace(0)), 1)
        self.assertEqual(len(mesh.CollapseEdge(0)), 2)

        mesh.AddFace([0, 1, 2])
        self.assertEqual(mesh.NumVertices(), 3)
        mesh.DisconnectFace(0)
        self.assertEqual(mesh.NumVertices(), 6)
        for _ in range(3):
            mesh.AddFace([0, 1, 2])
        self.assertEqual(mesh.NumVertices(), 6)
        mesh.DisconnectFacesUntilManifoldEdges()
        self.assertEqual(mesh.NumVertices(), 15)
        mesh.AddFace([0, mesh.VPosition(1), mesh.VPosition(2)])
        self.assertEqual(mesh.NumVertices(), 17)
        mesh.DisconnectFacesUntilManifoldVertices()
        self.assertEqual(mesh.NumVertices(), 18)
        mesh.AddFace([0, 1, 2])
        self.assertEqual(mesh.NumVertices(), 18)
        mesh.DisconnectFacesUntilManifold()
        self.assertEqual(mesh.NumVertices(), 21)
        for _ in range(mesh.NumFaces()):
            mesh.RemoveFace(0)
        self.assertEqual(mesh.NumVertices(), 0)
        self.assertEqual(mesh.NumFaces(), 0)
        self.assertEqual(mesh.NumHalfedges(), 0)


if __name__ == "__main__":
    sys.argv = sys.argv[:1]
    unittest.main()
