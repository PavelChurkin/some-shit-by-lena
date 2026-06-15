#!/usr/bin/env python3
"""Generate a self-contained static HTML visualization of the merge counts.

Reuses cube_merge_search.representations so the embedded numbers cannot drift
from the report. Run:  python3 "Теория Лены/generate_visualization.py"
Output: visualization.html (no external dependencies, opens offline).
"""
import json
import os

from cube_merge_search import representations

ZMAX = 40
KS = list(range(3, 9))

data = {
    "zmax": ZMAX,
    "ks": KS,
    "counts": {k: [len(representations(z, k)) for z in range(2, ZMAX + 1)]
               for k in KS},
    "zs": list(range(2, ZMAX + 1)),
}

HTML = """<!DOCTYPE html>
<html lang="ru">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Теория Лены — рост числа вариаций слияния</title>
<style>
  body { font-family: system-ui, sans-serif; margin: 2rem; background:#0f1115; color:#e6e6e6; }
  h1 { font-size: 1.3rem; }
  p  { max-width: 60ch; color:#b8b8b8; }
  #chart { background:#161922; border:1px solid #2a2f3a; border-radius:8px; }
  .legend span { display:inline-block; margin-right:1rem; font-size:.9rem; }
  .sw { display:inline-block; width:12px; height:12px; border-radius:2px; vertical-align:middle; margin-right:4px; }
  code { color:#9ad; }
</style>
</head>
<body>
<h1>Рост числа вариаций слияния: представления z³ суммой k кубов</h1>
<p>Чем больше масштаб z и чем больше число слагаемых k, тем больше вариантов
слияния. Канал <code>k = 8</code> — самый ёмкий. Данные вычислены скриптом
<code>cube_merge_search.py</code> и встроены в страницу.</p>
<div class="legend" id="legend"></div>
<canvas id="chart" width="900" height="480"></canvas>
<script>
const DATA = __DATA__;
const COLORS = {3:"#ff6b6b",4:"#ffa94d",5:"#ffd43b",6:"#69db7c",7:"#4dabf7",8:"#b197fc"};
const cv = document.getElementById("chart"), ctx = cv.getContext("2d");
const W = cv.width, H = cv.height, P = 50;
const zs = DATA.zs, maxC = Math.max(...DATA.ks.flatMap(k => DATA.counts[k]));
const x = z => P + (z - zs[0]) / (zs[zs.length-1] - zs[0]) * (W - 2*P);
const y = c => H - P - (c / maxC) * (H - 2*P);
// axes
ctx.strokeStyle = "#3a3f4a"; ctx.fillStyle = "#888"; ctx.font = "12px sans-serif";
ctx.beginPath(); ctx.moveTo(P,P); ctx.lineTo(P,H-P); ctx.lineTo(W-P,H-P); ctx.stroke();
ctx.textAlign = "right";
for (let i=0;i<=5;i++){ const c=Math.round(maxC*i/5); const yy=y(c);
  ctx.strokeStyle="#222630"; ctx.beginPath(); ctx.moveTo(P,yy); ctx.lineTo(W-P,yy); ctx.stroke();
  ctx.fillText(c, P-8, yy+4); }
ctx.textAlign = "left";
for (let z=zs[0]; z<=zs[zs.length-1]; z+=4){ ctx.fillText(z, x(z)-6, H-P+18); }
ctx.fillText("z (ребро итогового куба)", W/2-60, H-12);
ctx.save(); ctx.translate(14,H/2); ctx.rotate(-Math.PI/2); ctx.fillText("число вариантов слияния",0,0); ctx.restore();
// lines
for (const k of DATA.ks){ const cnts = DATA.counts[k];
  ctx.strokeStyle = COLORS[k]; ctx.lineWidth = 2; ctx.beginPath();
  zs.forEach((z,i)=>{ const px=x(z), py=y(cnts[i]); i?ctx.lineTo(px,py):ctx.moveTo(px,py); });
  ctx.stroke(); }
// legend
const leg = document.getElementById("legend");
DATA.ks.forEach(k => { const s=document.createElement("span");
  s.innerHTML = `<span class="sw" style="background:${COLORS[k]}"></span>k = ${k}`; leg.appendChild(s); });
</script>
</body>
</html>
"""

if __name__ == "__main__":
    out = os.path.join(os.path.dirname(__file__), "visualization.html")
    html = HTML.replace("__DATA__", json.dumps(data))
    with open(out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"wrote {out}")
