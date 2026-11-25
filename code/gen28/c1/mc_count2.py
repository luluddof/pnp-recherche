import sys
from multicolor_count import run
for h,m in [(500,25),(500,100),(500,500)]:
    run(h,m,sys.stdout)
