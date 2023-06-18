func orient2d*(px,  py,  qx,  qy,  rx,  ry: float): float =
  return (qy - py) * (rx - qx) - (qx - px) * (ry - qy)

