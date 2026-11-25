"""Modify only suffix storage allocation in official KMC 3.2.4 API.
Lookup, signatures, prefix decoding and counter parsing remain unchanged.
"""
from pathlib import Path
import shutil,hashlib,json,difflib
O=Path(__file__).resolve().parent;V=O/'vendor_kmc_api';B=O/'mmap_api';B.mkdir(exist_ok=True)
for p in V.glob('*'):
    if p.suffix in {'.h','.cpp'}:shutil.copy2(p,B/p.name)
h=(B/'kmc_file.h').read_text()
h=h.replace('protected:\n\tenum open_mode','protected:\n\tvoid* suffix_mapping = nullptr;\n\tsize_t suffix_mapping_size = 0;\n\tvoid ReleaseSuffix();\n\tenum open_mode',1)
assert 'suffix_mapping = nullptr' in h
(B/'kmc_file.h').write_text(h)
p=B/'kmc_file.cpp';s=p.read_text()
s=s.replace('#include <tuple>','#include <tuple>\n#include <sys/mman.h>\n#include <unistd.h>\n#include <stdexcept>',1)
old='\tsufix_file_buf = new uchar[size];\n\tresult = fread(sufix_file_buf, 1, size, file_suf);\n\tif (result != size)\n\t\treturn false;'
new='''\tsuffix_mapping_size = size + 8;
\tsuffix_mapping = mmap(nullptr, suffix_mapping_size, PROT_READ, MAP_PRIVATE, fileno(file_suf), 0);
\tif (suffix_mapping == MAP_FAILED) { suffix_mapping = nullptr; return false; }
\tmadvise(suffix_mapping, suffix_mapping_size, MADV_RANDOM);
\tsufix_file_buf = static_cast<uchar*>(suffix_mapping) + 4;'''
assert s.count(old)==1;s=s.replace(old,new)
old='\tif (sufix_file_buf)\n\t\tdelete[] sufix_file_buf;'
assert s.count(old)==1;s=s.replace(old,'\tReleaseSuffix();')
old='\t\tdelete [] sufix_file_buf;\n\t\tsufix_file_buf = NULL;'
assert s.count(old)==1;s=s.replace(old,'\t\tReleaseSuffix();')
s+='''
void CKMCFile::ReleaseSuffix() {
    if (suffix_mapping) {
        munmap(suffix_mapping, suffix_mapping_size);
        suffix_mapping = nullptr;
        suffix_mapping_size = 0;
    } else if (sufix_file_buf) { delete[] sufix_file_buf; }
    sufix_file_buf = nullptr;
}
'''
p.write_text(s)
patch=''
for name in ['kmc_file.h','kmc_file.cpp']:
    patch+=''.join(difflib.unified_diff((V/name).read_text().splitlines(True),(B/name).read_text().splitlines(True),fromfile='official/'+name,tofile='mmap/'+name))
(O/'mmap_storage.patch').write_text(patch)
(O/'vendor_manifest.json').write_text(json.dumps({p.name:dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest(),source='https://github.com/refresh-bio/KMC/tree/v3.2.4/kmc_api') for p in V.iterdir() if p.is_file()},indent=2)+'\n')
print('Patched suffix storage only; full diff saved')
