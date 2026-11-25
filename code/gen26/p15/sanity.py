import json,sys
sys.path.insert(0,'.')
import importlib.util
spec=importlib.util.spec_from_file_location("mv","myverify.py")
src=open("myverify.py").read().split("claims=")[0]; ns={}; exec(src,ns)
# C=0 -> rank N^2 ; perturb one entry of a valid file -> rank should rise
d=json.load(open("/home/user/pnp-recherche/code/gen26/p11/exact_6_8_6.json"))
z={"N":6,"E":[[[0]*6 for _ in range(6)]],"C":[[[0]*6 for _ in range(6)]]}
json.dump(z,open("z.json","w")); print("C=0:",ns["check"]("z.json",0,1))
d["E"][0][0][0]+=1; json.dump(d,open("pert.json","w")); print("perturbed:",ns["check"]("pert.json",8,6))
