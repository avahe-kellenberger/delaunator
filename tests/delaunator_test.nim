import delaunator, nimtest

import fixtures/issue13

proc validate(coords: seq[float])

describe "Delaunator":

  it "produces correct triangulation":
    validate(data)

proc validate(coords: seq[float]) =
  let d = newDelaunator(coords)
  # Validate halfEdges
  echo d.halfEdges.len
  for i in 0 ..< d.halfEdges.len:
    # TODO: Seems like the last half edge is wrong!
    if d.halfEdges[i] != -1 and d.halfEdges[d.halfEdges[i]] != i:
      echo $i
      echo $d.halfEdges[i]
      echo $d.halfEdges[d.halfEdges[i]]
      break

