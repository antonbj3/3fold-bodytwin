from pathlib import Path
import requests,json,hashlib,datetime,concurrent.futures,subprocess,time
P=Path(__file__).resolve().parent/'literature/r4';P.mkdir(exist_ok=True)
SOURCES={
'barnett2016':'https://citeseerx.ist.psu.edu/document?doi=0f2487feb9d45d67e4e258e3f900c0d761a28fbc&repid=rep1&type=pdf',
'owen2022':'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC8876719/unicode',
'barnett_thesis_search':'https://etda.libraries.psu.edu/catalog?q=Barnett+Andrew&search_field=all_fields',
'irwin2021_openalex':'https://api.openalex.org/works/https://doi.org/10.1016/j.jmbbm.2021.104660',
'barnett2016_openalex':'https://api.openalex.org/works/https://doi.org/10.1115/1.4030374'}
def fetch(item):
 name,url=item;start=time.perf_counter();target=P/(name+'.response');meta={'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 if target.exists():return {'name':name,'already_preserved':True}
 try:
  r=requests.get(url,timeout=35);target.open('xb').write(r.content);meta.update(status=r.status_code,final_url=r.url,bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest(),is_pdf=r.content.startswith(b'%PDF'))
  if meta['is_pdf']:
   subprocess.run(['pdftotext','-layout',str(target),str(P/(name+'.txt'))],check=True)
 except Exception as e:meta['error']=str(e)
 meta['wall_s']=time.perf_counter()-start;(P/(name+'.acquisition.json')).open('x').write(json.dumps(meta,indent=2)+'\n');return {'name':name,**meta}
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
  for v in ex.map(fetch,SOURCES.items()):print(json.dumps(v))
