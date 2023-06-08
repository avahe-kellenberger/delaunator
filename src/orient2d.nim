func orient2d*(px,  py,  qx,  qy,  rx,  ry: float): bool =
  return (qy - py) * (rx - qx) - (qx - px) * (ry - qy) < 0;

