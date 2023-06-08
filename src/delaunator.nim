import std/math
import orient2d as orient2dModule

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

    cx: float
    cy: float
    hullStart: int
    trianglesLen: int

proc update(this: Delaunator)
func hashKey(this: Delaunator, x, y: int): int
proc addTriangle(this: Delaunator, i0, i1, i2, a, b, c: int): int
proc link(this: Delaunator, a, b: int)

func dist(ax, ay, bx, by: float): float
func circumradius(ax, ay, bx, by, cx, cy: float): float
func circumcenter(ax,  ay,  bx,  by,  cx,  cy: float): Vector
func quicksort(ids: var seq[uint32], dists: var seq[float64], left, right: int)
func swap[T: SomeNumber](arr: var seq[T], i, j: int)
proc fill[T](arr: seq[T], val: T)

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
  let n = this.coords.len shl 1

  var
    minX = Inf
    minY = Inf
    maxX = NegInf
    maxY = NegInf

  for i in 0 ..< n:
    let x = this.coords[2 * i]
    let y = this.coords[2 * i + 1]
    minX = min(x, minX)
    minY = min(y, minY)
    maxX = max(x, maxX)
    maxY = max(y, maxY)
    this.ids[i] = i

  let cx = (minX + maxX) / 2
  let cy = (minY + maxY) / 2

  var
    i0: int
    i1: int
    i2: int

  # Using a block to prevent scope pollution (minDist variable)
  block:
    # Pick a seed point close to the center
    var minDist = Inf
    for i in 0 ..< n:
      let d = dist(cx, cy, this.coords[2 * i], this.coords[2 * i + 1])
      if d < minDist:
        i0 = i
        minDist = d

  let i0x = this.coords[2 * i0]
  let i0y = this.coords[2 * i0 + 1]
  
  block:
    # Find the point closest to the seed
    var minDist = Inf
    for i in 0 ..< n:
      if i == i0:
        continue
      
      let d = dist(i0x, i0y, this.coords[2 * i], this.coords[2 * i + 1])
      if d < minDist and d > 0:
        i1 = i
        minDist = d

  let i1x = this.coords[2 * i1]
  let i1y = this.coords[2 * i1 + 1]

  # Find the third point which forms the smallest circumcircle with the first two
  var minRadius = Inf
  for i in 0 ..< n:
    if i == i0 or i == i1:
      continue

    let r = circumradius(i0x, i0y, i1x, i1y, this.coords[2 * i], this.coords[2 * i + 1])
    if r < minRadius:
      i2 = i
      minRadius = r

  let i2x = this.coords[2 * i2]
  let i2y = this.coords[2 * i2 + 1]

  if minRadius == Inf:
    # Order collinear points by dx (or dy if all x are identical)
    # and return the list as a hull
    for i in 0 ..< n:
      let x = this.coords[2 * i] - this.coords[0]
      if x != 0:
        this.dists[i] = x
      else:
        this.dists[i] = this.coords[2 * i + 1] - this.coords[1]

    quicksort(this.ids, this.dists, 0, n - 1)

    let hull = newSeqOfCap[uint32](n)
    var
      j = 0
      d0 = NegInf
    for i in 0 ..< n:
      let id = this.ids[i]
      let d = this.dists[id]
      if d > d0:
        hull[j] = id
        d0 = d
        inc j

    this.hull = hull[0 .. j - 1]
    this.triangles = newSeq[uint32]()
    this.halfEdges = newSeq[uint32]()
    return

  # Swap the order of the seed points for counter-clockwise orientation
  if orient2d(i0x, i0y, i1x, i1y, i2x, i2y):
    let
      i = i1
      x = i1x
      y = i1y

    i1 = i2
    i1x = i2x
    i1y = i2y
    i2 = i
    i2x = x
    i2y = y

  let center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y)
  this.cx = center.x
  this.cy = center.y

  for i in 0 ..< n:
    this.dists[i] = dist(this.coords[2 * i], this.coords[2 * i + 1], center.x, center.y)

  # Sort the points by distance from the seed triangle circumcenter
  quicksort(this.ids, this.dists, 0, n - 1)

  # Set up the seed triangle as the starting hull
  this.hullStart = i0
  let hullSize = 3
  
  this.hullNext[i0] = i1
  this.hullNext[i1] = i2
  this.hullNext[i2] = i0

  this.hullPrev[i2] = i1
  this.hullPrev[i0] = i2
  this.hullPrev[i1] = i0

  this.hullTri[i0] = 0
  this.hullTri[i1] = 1
  this.hullTri[i2] = 2

  this.hullHash.fill(-1)
  this.hullHash[this.hashKey(i0x, i0y)] = i0
  this.hullHash[this.hashKey(i1x, i1y)] = i1
  this.hullHash[this.hashKey(i2x, i2y)] = i2

  this.trianglesLen = 0
  this.addTriangle(i0, i1, i2, -1, -1, -1)

  var xp, yp: int
  for k in 0 ..< this.ids.len:
    let
      i = this.ids[k]
      x = this.coords[2 * i]
      y = this.coords[2 * i + 1]

    if k > 0 and abs(x - xp) <= EPSILON and abs(y - yp) <= EPSILON:
      continue

proc link(this: Delaunator, a, b: int) =
  this.halfEdges[a] = b
  if b != -1:
    this.halfEdges[b] = a

proc addTriangle(this: Delaunator, i0, i1, i2, a, b, c: int): int =
  ## Add a new triangle given vertex indices and adjacent half-edge ids
  result = this.trianglesLen

  this.triangles[result] = i0
  this.triangles[result + 1] = i1
  this.triangles[result + 2] = i2

  this.link(result, a)
  this.link(result + 1, b)
  this.link(result + 2, c)

  this.trianglesLen += 3

func dist(ax, ay, bx, by: float): float =
  let
    dx = ax - bx
    dy = ay - by
  return dx * dx + dy * dy

func circumradius(ax, ay, bx, by, cx, cy: float): float =
  let dx = bx - ax
  let dy = by - ay
  let ex = cx - ax
  let ey = cy - ay

  let bl = dx * dx + dy * dy
  let cl = ex * ex + ey * ey
  let d = 0.5 / (dx * ey - dy * ex)

  let x = (ey * bl - dy * cl) * d
  let y = (dx * cl - ex * bl) * d

  return x * x + y * y

func circumcenter(ax,  ay,  bx,  by,  cx,  cy: float): Vector =
  let dx = bx - ax
  let dy = by - ay
  let ex = cx - ax
  let ey = cy - ay
  let bl = dx * dx + dy * dy
  let cl = ex * ex + ey * ey
  let d = 0.5 / (dx * ey - dy * ex)
  let x = ax + (ey * bl - dy * cl) * d
  let y = ay + (dx * cl - ex * bl) * d
  return (x, y)

template doWhile(a, b: untyped): untyped =
  b
  while a:
    b

func quicksort(ids: var seq[uint32], dists: var seq[float64], left, right: int) =
  if (right - left) <= 20:
    for i in (left + 1) .. right:
      let temp = ids[i]
      let tempDist = dists[temp]
      var j = i - 1
      while j >= left and dists[ids[j]] > tempDist:
        ids[j + 1] = ids[j]
        dec j
      ids[j + 1] = temp
  else:
    let median = (left + right) shr 1
    var i = left + 1
    var j = right
    swap(ids, median, i)
    if dists[ids[left]] > dists[ids[right]]:
      swap(ids, left, right)
    if dists[ids[i]] > dists[ids[right]]:
      swap(ids, i, right)
    if dists[ids[left]] > dists[ids[i]]:
      swap(ids, left, i)

    let temp = ids[i]
    let tempDist = dists[temp]
    while true:
      doWhile (dists[ids[i]] < tempDist): inc i
      doWhile (dists[ids[j]] > tempDist): dec j
      if j < i:
        break
      swap(ids, i, j)
    
    ids[left + 1] = ids[j]
    ids[j] = temp

    if (right - i + 1) >= (j - left):
      quicksort(ids, dists, i, right)
      quicksort(ids, dists, left, j - 1)
    else:
      quicksort(ids, dists, left, j - 1)
      quicksort(ids, dists, i, right)

func swap[T: SomeNumber](arr: var seq[T], i, j: int) =
  let tmp = arr[i]
  arr[i] = arr[j]
  arr[j] = tmp

proc fill[T](arr: seq[T], val: T) =
  for i in 0 ..< arr:
    arr[i] = val

func pseudoAngle(dx, dy: float): float =
  ## Monotonically increases with real angle, but doesn't need expensive trigonometry
  let p = dx / (abs(dx) + abs(dy))
  result =
    if dy > 0:
      3 - p
    else:
      1 + p

  result /= 4

func hashKey(this: Delaunator, x, y: int): int =
  return int floor(pseudoAngle(x - this.cx, y - this.cy) * this.hashSize) mod this.hashSize

