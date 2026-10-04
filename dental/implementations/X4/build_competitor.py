"""Compile the pinned, unmodified OsteoOpt geometric simplifier.

Minimal vector adapter implements exactly the seven Euclidean operations used
by this file. It does not implement or replace ArtiSynth/union optimization.
Upstream PolyForm noncommercial license and NOTICE are saved in sources/.
"""
import shutil, subprocess
from pathlib import Path
from freeze import ROOT, write, sha

def build():
    base = ROOT / 'comparison_java'
    src = base / 'artisynth/VSP/reconstruction'
    src.mkdir(parents=True, exist_ok=True)
    original = ROOT / 'sources/osteoopt_PolylineSimplifier.java'
    target = src / 'PolylineSimplifier.java'
    target.write_bytes(original.read_bytes())
    assert sha(original) == sha(target)
    matrix = base / 'maspack/matrix'
    matrix.mkdir(parents=True, exist_ok=True)
    (matrix / 'Vector3d.java').write_text('package maspack.matrix;\npublic class Vector3d {\n public double x,y,z;\n public Vector3d() {}\n public Vector3d(Vector3d v) {set(v);}\n public Vector3d(double a,double b,double c) {x=a;y=b;z=c;}\n public void set(Vector3d v) {x=v.x;y=v.y;z=v.z;}\n public void sub(Vector3d a,Vector3d b) {x=a.x-b.x;y=a.y-b.y;z=a.z-b.z;}\n public double normSquared() {return x*x+y*y+z*z;}\n public double dot(Vector3d b) {return x*b.x+y*b.y+z*b.z;}\n public void scaledAdd(double s,Vector3d b) {x+=s*b.x;y+=s*b.y;z+=s*b.z;}\n public double distanceSquared(Vector3d b) {double a=x-b.x,c=y-b.y,d=z-b.z;return a*a+c*c+d*d;}\n}')
    (matrix / 'Point3d.java').write_text('package maspack.matrix;\npublic class Point3d extends Vector3d {\n public Point3d(Vector3d p) {super(p);}\n public Point3d(double x,double y,double z) {super(x,y,z);}\n public double distance(Point3d p) {return Math.sqrt(distanceSquared(p));}\n}')
    (base / 'Compare.java').write_text('import java.io.*;\nimport java.util.*;\nimport maspack.matrix.Point3d;\nimport artisynth.VSP.reconstruction.PolylineSimplifier;\npublic class Compare {\n public static void main(String[] args) throws Exception {\n  Scanner s=new Scanner(System.in);s.useLocale(Locale.US);\n  int n=s.nextInt(),k=s.nextInt();double min=s.nextDouble();\n  ArrayList<Point3d> curve=new ArrayList<>();\n  for(int i=0;i<n;i++)curve.add(new Point3d(s.nextDouble(),s.nextDouble(),s.nextDouble()));\n  ArrayList<Point3d> nodes=PolylineSimplifier.bisectSimplifyDouglasPeucker(curve,min,k);\n  System.out.println("NODES "+nodes.size());\n  for(Point3d p:nodes)System.out.printf(Locale.US,"%.17g %.17g %.17g%n",p.x,p.y,p.z);\n }\n}')
    files = [str(p) for p in base.rglob('*.java')]
    result = subprocess.run(['javac', *files], capture_output=True, text=True)
    write(ROOT / 'COMPETITOR_BUILD.json', dict(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr, upstream_class_sha256=sha(original), copied_class_sha256=sha(target), javac=shutil.which('javac'), scope='Unmodified upstream geometric simplifier; Euclidean-vector adapter; full OsteoOpt++ BO/ArtiSynth NOT_RUN', full_system_missing=['MATLAB executable', 'ArtiSynth runtime', 'matched donor geometry', 'CT parameters and loads', 'patient-specific meshes upstream not released'], source_manifest=str(ROOT / 'sources/osteoopt_manifest.json')))
    if result.returncode:
        raise RuntimeError(result.stderr)
if __name__ == '__main__':
    build()
