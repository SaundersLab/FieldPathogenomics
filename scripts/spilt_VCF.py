import sys,os, io
import subprocess

base = sys.argv[1]
N_split = int(sys.argv[2])

genome_length = 63e6
per_splt = int(genome_length/N_split)

print("Using {0} records per file".format(per_splt))

def rotate_output(base, i):
    proc = subprocess.Popen('/tgac/software/production/tabix/0.2.6/x86_64/bin/bgzip -c > {0}'.format(base+'_'+str(i)+".vcf.gz"), shell=True, stdin=subprocess.PIPE, bufsize=1, universal_newlines=True)
    return proc

if __name__ == '__main__':
    got_header = False
    header = ''
    record_count, split_count = 0,0
    buf = io.StringIO()
    gzip_proc = []

    for line in sys.stdin:
        
        # Catch the header
        if got_header is False:
            if line[0] == '#':
                header+=line
                continue
            elif header is '':
                raise Exception("No header on input VCF stream")
            else:
                got_header = True
                gzip_proc.append(rotate_output(base, split_count))
                gzip_proc[split_count].stdin.write(header)
                
        gzip_proc[split_count].stdin.write(line)
        record_count+=1
        
        if record_count > per_splt and split_count < N_split:
            # If we've filled this split and still have some left
            # Flush to this gzip instance and send eof so it starts compressing
            gzip_proc[split_count].stdin.flush()
            gzip_proc[split_count].stdin.close()
            
            print("Compressing file {0}".format(base+'_'+str(split_count)+".vcf.gz"))
            # Set up the next split
            split_count+=1
            record_count = 0
            gzip_proc.append(rotate_output(base, split_count))
            gzip_proc[split_count].stdin.write(header)
    
    # Now tidy up the final split
    gzip_proc[split_count].stdin.flush()
    gzip_proc[split_count].stdin.close()
    
    #Wait for a gzip procs to finsh sucessfully
    for i,proc in enumerate(gzip_proc):
        proc.wait()
        if proc.returncode != 0:
            raise Exception("BZGIP process for file {0} failed".format(base+'_'+str(i)+".vcf.gz"))
            
    