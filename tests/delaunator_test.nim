import delaunator, nimtest

import fixtures/issue13

proc validate(coords: seq[float])

describe "Delaunator":

  it "produces correct triangulation":
    validate(data)

proc validate(coords: seq[float]) =
  let d = newDelaunator(coords)
  # Validate halfEdges
  for i in 0 ..< d.halfEdges.len:
    # if d.halfEdges[i] != -1 and d.halfEdges[d.halfEdges[i]] != i:
    if d.halfEdges[i] != -1 and d.halfEdges[d.halfEdges[i]] != i:
      echo "Wrong: "
      echo "i: " & $i
      echo "halfEdges[i]: " & $d.halfEdges[i]
      echo "halfEdges[halfEdges[i]]: " & $d.halfEdges[d.halfEdges[i]]
      echo ""


