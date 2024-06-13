# MIT License
#
# Copyright (c) 2024 TriLite https://github.com/MeshLite/TriLite
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

    def test_basic_properties(self):
        mesh = self.__class__.mesh
        self.assertEqual(mesh.NumHalfedges(), 6)
        self.assertEqual(mesh.NumVertices(), 4)
        self.assertEqual(mesh.NumFaces(), 2)

    def test_halfedges_vertices_faces(self):
        mesh = self.__class__.mesh
        self.assertEqual(len(mesh.Halfedges()), 6)
        self.assertEqual(len(mesh.Vertices()), 4)
        self.assertEqual(len(mesh.Faces()), 2)

    def test_positions(self):
        mesh = self.__class__.mesh
        self.assertEqual(len(mesh.Positions()), 4)

    def test_halfedge_navigation(self):
        mesh = self.__class__.mesh
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

    def test_geometric_properties(self):
        mesh = self.__class__.mesh
        self.assertTrue(
            np.array_equal(mesh.HGeometry(0), np.array([1.0, 0.0, 0.0]))
        )
        self.assertTrue(
            np.array_equal(mesh.HCentroid(1), np.array([1.0, 0.5, 0.0]))
        )
        self.assertEqual(len(mesh.HConnectionsAroundStart(0)), 2)
        self.assertEqual(len(mesh.HHalfedgesAroundHole(3)), 4)
        self.assertAlmostEqual(mesh.HLength(0), 1.0)

    def test_vertex_properties(self):
        mesh = self.__class__.mesh
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
        self.assertTrue(mesh.VIsBoundary(0))
        self.assertAlmostEqual(mesh.MedianEdgeLength(), 1.0)

    def test_face_properties(self):
        mesh = self.__class__.mesh
        self.assertEqual(mesh.FHalfedge(0), 0)
        self.assertEqual(len(mesh.FHalfedges(0)), 3)
        self.assertEqual(mesh.FNeighbors(0), [1])
        self.assertEqual(len(mesh.FVertices(0)), 3)
        self.assertEqual(len(mesh.FPositions(0)), 3)
        self.assertTrue(np.allclose(mesh.FNormal(0), [0, 0, 1]))
        self.assertAlmostEqual(mesh.FArea(0), 0.5)

    def test_edge_properties(self):
        mesh = self.__class__.mesh
        self.assertEqual(len(mesh.EdgeHalfedges(0)), 1)
        self.assertEqual(len(mesh.EdgeFaces(0)), 1)
        self.assertTrue(mesh.EdgeIsManifold(0))

    def test_boundary_halfedges(self):
        mesh = self.__class__.mesh
        self.assertGreaterEqual(len(mesh.BoundaryHalfedges()), 0)

    def test_centroid(self):
        mesh = self.__class__.mesh
        f_centroid_0 = mesh.FCentroid(0)
        expected_f_centroid_0 = np.array([2 / 3, 1 / 3, 0.0])
        self.assertTrue(np.allclose(f_centroid_0, expected_f_centroid_0))
        centroid = mesh.Centroid()
        expected_centroid = np.array([0.5, 0.5, 0.0])
        self.assertTrue(np.allclose(centroid, expected_centroid))

    def test_bounding_box(self):
        mesh = self.__class__.mesh
        f_bbox_0 = mesh.FBoundingBox(0)
        expected_f_bbox_0 = (
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 1.0, 0.0]),
        )
        self.assertTrue(np.allclose(f_bbox_0[0], expected_f_bbox_0[0]))
        self.assertTrue(np.allclose(f_bbox_0[1], expected_f_bbox_0[1]))
        bbox = mesh.BoundingBox()
        expected_bbox = (np.array([0.0, 0.0, 0.0]), np.array([1.0, 1.0, 0.0]))
        self.assertTrue(np.allclose(bbox[0], expected_bbox[0]))
        self.assertTrue(np.allclose(bbox[1], expected_bbox[1]))

    def test_modify_mesh(self):
        mesh = TL.Trimesh(self.__class__.mesh)
        mesh.AddFace(
            [
                np.array([0.0, 1.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([1.0, 1.0, 0.0]),
            ]
        )
        self.assertEqual(len(mesh.RemoveFace(0)), 1)
        self.assertEqual(
            mesh.CollapseEdge(0, mesh.HCentroid(0)), ([0], [1, 1, 1])
        )
        self.assertEqual(mesh.RemoveFace(0), [0, 1, 0])
        mesh.AddFace(
            [
                np.array([0.0, 1.0, 0.0]),
                np.array([1.0, 0.0, 0.0]),
                np.array([1.0, 1.0, 0.0]),
            ]
        )
        mesh.AddFace([mesh.HEnd(0), mesh.HStart(0), np.array([1.0, 1.0, 1.0])])
        mesh.SplitEdge(0)
        self.assertEqual(mesh.NumFaces(), 4)
        mesh.CollapseEdge(0, mesh.HCentroid(0))
        self.assertEqual(mesh.NumFaces(), 2)
        mesh.FlipHalfedgeWithOpposite(0)
        self.assertEqual(mesh.HStart(0), 3)
        for _ in range(3):
            mesh.FlipHalfedgeWithOpposite(0)
        self.assertEqual(mesh.HStart(0), 0)
        mesh.AddFace([0, 1, 2])
        self.assertEqual(mesh.FNeighbors(0), [1])
        mesh.RemoveFace(1)
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

        v0 = np.array([0.0, 0.0, 0.0])
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([1.0, 1.0, 0.0])
        v3 = np.array([0.0, 1.0, 0.0])

        mesh.AddFace([v0, v1, v2])
        mesh.AddFace([2, v3, 0])

    def test_transformations(self):
        mesh = TL.Trimesh(self.__class__.mesh)
        copy = TL.Trimesh(mesh)
        rotation_matrix = np.array(
            [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
        )
        mesh.ApplyRotation(rotation_matrix)
        for v in range(mesh.NumVertices()):
            pos = mesh.VPosition(v)
            self.assertTrue(
                np.allclose(pos, rotation_matrix @ copy.VPosition(v))
            )

        copy = TL.Trimesh(mesh)
        translation_vector = np.array([1.0, 1.0, 1.0])
        mesh.ApplyTranslation(translation_vector)
        for v in range(mesh.NumVertices()):
            pos = mesh.VPosition(v)
            self.assertTrue(
                np.allclose(pos, copy.VPosition(v) + translation_vector)
            )

        copy = TL.Trimesh(mesh)
        scaling_factor = 2.0
        mesh.ApplyScaling(scaling_factor)
        for v in range(mesh.NumVertices()):
            pos = mesh.VPosition(v)
            self.assertTrue(np.allclose(pos, 2 * copy.VPosition(v)))

        copy = TL.Trimesh(mesh)
        mesh.ApplyRigidTransformation(rotation_matrix, translation_vector)
        for v in range(mesh.NumVertices()):
            pos = mesh.VPosition(v)
            self.assertTrue(
                np.allclose(
                    pos,
                    (rotation_matrix @ copy.VPosition(v)) + translation_vector,
                )
            )

        # Apply and test similarity transformation (rotation + scaling + translation)
        copy = TL.Trimesh(mesh)
        rotation_matrix = np.eye(3)
        translation_vector = np.array([1.0, 1.0, 1.0])
        scaling_factor = 2.0
        mesh.ApplySimilarityTransformation(
            rotation_matrix, translation_vector, scaling_factor
        )
        for v in range(mesh.NumVertices()):
            pos = mesh.VPosition(v)
            self.assertTrue(
                np.allclose(
                    pos,
                    scaling_factor * (rotation_matrix @ copy.VPosition(v))
                    + translation_vector,
                )
            )


if __name__ == "__main__":
    sys.argv = sys.argv[:1]
    unittest.main()
