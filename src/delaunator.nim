import std/math

const EPSILON = pow(2.0, -52.0)

var EDGE_STACK: array[512, uint32]

type
  Vector* = concept v
    v.x is float
    v.y is float

  Delaunator*[T: Vector] = object
    coords: seq[T]

    # Arrays that will store the triangulation graph
    triangles: seq[uint32]
    halfEdges: seq[int32]

    # Temporary arrays for tracking the edges of the advancing convex hull
    hashSize: int
    hullPrev: seq[uint32]
    hullNext: seq[uint32]
    hullTri: seq[uint32]
    hullHash: seq[int32]

    # Temporary arrays for sorting points
    ids: seq[uint32]
    dists: seq[float64]

proc update(this: Delaunator)

proc newDelaunator*[T: Vector](coords: seq[T]): Delaunator[T] =
  let n = coords.len shl 1
  
  result.coords = coords

  let maxTriangles = max(2 * n - 5, 0)
  result.triangles = newSeqOfCap[uint32](maxTriangles * 3)
  result.halfEdges = newSeqOfCap[int32](maxTriangles * 3)

  result.hashSize = ceil(sqrt(n))
  result.hullPrev = newSeqOfCap[uint32](n)
  result.hullNext = newSeqOfCap[uint32](n)
  result.hullTri = newSeqOfCap[uint32](n)
  result.hullHash = newSeqOfCap[int32](result.hashSize)

  result.ids = newSeqOfCap[uint32](n)
  result.dists = newSeqOfCap[float64](n)

  result.update()

proc update(this: Delaunator) =
  # TODO
  discard

