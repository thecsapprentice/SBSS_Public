#include "ARCHIVE.h"

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>

#include "libtar.h"
#include <errno.h>
#include <cstring>
#include <cstdio>
#include <fcntl.h>
#include <sstream>
#include <fstream>

using namespace PhysBAM;

namespace ARCHIVE_WORKSPACE{

    typedef struct{std::string filename; ARRAY<char> data; int last_byte; bool open;} memory_file;
    
    std::ostream& operator<<(std::ostream& out, const memory_file& mf)
    {out << mf.filename << ": "<< mf.data.m << "B"; return out;}
    
    typedef ARRAY<memory_file> file_registry;
    file_registry files;
        
    int open_memoryfile(const char * filename, int oflags, ...)
    {
        //LOG::SCOPE scope("open_memoryfile");
        std::string s_filename(filename);
        //LOG::cout << "Attempting to open file " << s_filename << std::endl;
        
        int i;
        for(i=1;i<=files.m;i++){
            if(s_filename == files(i).filename){
                if(!files(i).open){
                    files(i).open=true;
                    return i;
                }
                else
                    return -1;
            }
        }     
        
        memory_file newFile;
        newFile.filename = s_filename;
        newFile.last_byte=1;
        newFile.open=true;
        
        return files.Append(newFile);
    }
    
    int close_memoryfile(int fidesc)
    {
        //LOG::SCOPE scope("close_memoryfile");
        
        if(fidesc == 0 || fidesc > files.m)
            return -1;
        
        if(!files(fidesc).open)
            return -1;
        
        files(fidesc).last_byte=1;
        files(fidesc).open=false;
        return 0;
    }
    
    ssize_t read_memoryfile(int fidesc, void * buf, size_t nbytes)
    {
        //LOG::SCOPE scope("read_memoryfile");
        
        if(fidesc == 0 || fidesc > files.m)
            return -1;
        if(!files(fidesc).open)
            return -1;
        
        memory_file& file = files(fidesc);
        
        char* cbuf = (char*)buf;
        
        int byte;
        for(byte=0;byte<nbytes && (byte+file.last_byte) <= file.data.m;byte++){
            cbuf[byte] = file.data(byte+file.last_byte);
        }
        file.last_byte+=byte;
        
        return byte;
    }
    
    ssize_t write_memoryfile(int fidesc, const void * buf, size_t nbytes)
    {
        //LOG::SCOPE scope("write_memoryfile");
        
        if(fidesc == 0 || fidesc > files.m)
            return -1;
        if(!files(fidesc).open)
            return -1;
        
        memory_file& file = files(fidesc);
        
        char* cbuf = (char*)buf;
        
        int byte;
        file.data.Resize(file.data.m + nbytes);
        //LOG::cout << "data buffer: " << file.data.m << std::endl;
        for(byte=0;byte<nbytes;byte++){
            //LOG::cout << "Writing byte : " << byte+file.last_byte << std::endl;
            file.data(byte+file.last_byte) = cbuf[byte];
        }
        file.last_byte+=byte;
        
        return byte;
    }
    
/* add file contents to a tarchive */
    int
    tar_append_memfile(TAR *t, const char *realname)
    {
        char block[T_BLOCKSIZE];
        int filefd;
        int i, j;
        size_t size;
        
        filefd = open_memoryfile(realname, O_RDONLY);
        if (filefd == -1)
            return -1;
        
        size = th_get_size(t);
        for (i = size; i > T_BLOCKSIZE; i -= T_BLOCKSIZE)
            {
                j = read_memoryfile(filefd, &block, T_BLOCKSIZE);
                if (j != T_BLOCKSIZE)
                    {
                        if (j != -1)
                            errno = EINVAL;
                        return -1;
                    }
                if (tar_block_write(t, &block) == -1)
                    return -1;
            }
        
        if (i > 0)
            {
                j = read_memoryfile(filefd, &block, i);
                if (j == -1)
                    return -1;
                memset(&(block[i]), 0, T_BLOCKSIZE - i);
                if (tar_block_write(t, &block) == -1)
			return -1;
            }
        
        close_memoryfile(filefd);
        
	return 0;
    }

    tartype_t memory_tar_type;
   
}


//***************************************************************
//
//        Public Functions 
//
//***************************************************************

void* ARCHIVE::archive_ptr=NULL;

void ARCHIVE::CreateArchive(const std::string archive_name)
{
    //LOG::SCOPE scope("ARCHIVE::CreateArchive");

    if( archive_ptr )
        PHYSBAM_FATAL_ERROR("Can't work on new archive. Archive already open!");

    ARCHIVE_WORKSPACE::memory_tar_type.openfunc = &ARCHIVE_WORKSPACE::open_memoryfile;
    ARCHIVE_WORKSPACE::memory_tar_type.closefunc = &ARCHIVE_WORKSPACE::close_memoryfile;
    ARCHIVE_WORKSPACE::memory_tar_type.readfunc = &ARCHIVE_WORKSPACE::read_memoryfile;
    ARCHIVE_WORKSPACE::memory_tar_type.writefunc = &ARCHIVE_WORKSPACE::write_memoryfile;

    char aname[512];
    strcpy( aname, archive_name.c_str() );

    if( tar_open((TAR**)(&archive_ptr), aname,
                 &ARCHIVE_WORKSPACE::memory_tar_type, O_WRONLY, 0, 0) ){
        PHYSBAM_FATAL_ERROR("Error opening In-Memory archive.");
    }
   
}

void ARCHIVE::OpenArchive(std::string archive_name)
{
    //LOG::SCOPE scope("ARCHIVE::OpenArchive");

    if( archive_ptr )
        PHYSBAM_FATAL_ERROR("Can't work on new archive. Archive already open!");
    
    char aname[512];
    strcpy( aname, archive_name.c_str() );

    if( tar_open((TAR**)(&archive_ptr), aname,
                 0, O_RDONLY, 0, 0) ){
        PHYSBAM_FATAL_ERROR("Error opening File archive.");
    }


}

void ARCHIVE::WriteDataBuffer(std::string identifier, const DATA_BUFFER& data)
{
    //LOG::SCOPE scope("ARCHIVE::WriteDataBuffer");

    if(!archive_ptr)
        PHYSBAM_FATAL_ERROR("Can't write archive. Archive not open!");

    TAR* tar_ptr = (TAR*)(archive_ptr);

    int mem_file = ARCHIVE_WORKSPACE::open_memoryfile(identifier.c_str(), O_WRONLY);
    ARCHIVE_WORKSPACE::write_memoryfile( mem_file, data.Get_Array_Pointer() , data.m );
    ARCHIVE_WORKSPACE::close_memoryfile( mem_file );

    char aname[512];
    strcpy( aname, identifier.c_str() );

    th_set_type( tar_ptr, S_IFREG);
    th_set_path( tar_ptr, aname );
    th_set_link( tar_ptr, aname );
    th_set_user( tar_ptr, 0 );
    th_set_group( tar_ptr, 0 );
    th_set_mode( tar_ptr, S_IRWXU );
    th_set_size( tar_ptr, data.m);
    th_write(tar_ptr);
    ARCHIVE_WORKSPACE::tar_append_memfile(tar_ptr, identifier.c_str() );
    th_finish( tar_ptr );
}


void ARCHIVE::ReadDataBuffer(std::string identifier, DATA_BUFFER& data)
{
    //LOG::SCOPE scope("ARCHIVE::ReadDataBuffer");

    if(!archive_ptr)
        PHYSBAM_FATAL_ERROR("Can't read archive. Archive not open!");

    char buf[T_BLOCKSIZE];
    TAR* tar_ptr = (TAR*)(archive_ptr);

    //rewind((FILE*)(tar_ptr->fd));
    th_read(tar_ptr);
    do{
        if( strcmp(identifier.c_str(), th_get_pathname(tar_ptr))==0 ){
            int mem_file = ARCHIVE_WORKSPACE::open_memoryfile(identifier.c_str(), O_WRONLY);
            int size = th_get_size(tar_ptr);
            /* extract the file */
            for (int i = size; i > 0; i -= T_BLOCKSIZE)
                {
                    int k = tar_block_read(tar_ptr, buf);
                    if (k != T_BLOCKSIZE)
                        {
                            if (k != -1)
                                errno = EINVAL;
                            PHYSBAM_FATAL_ERROR("Can't read archive. Read Error Occured!");
                        }

                    /* write block to output file */
                    if (ARCHIVE_WORKSPACE::write_memoryfile(mem_file, buf,
                                                            ((i > T_BLOCKSIZE) ? T_BLOCKSIZE : i)) == -1)
                        PHYSBAM_FATAL_ERROR("Can't read archive. Write to Memory File Failed!");

                }

            ARCHIVE_WORKSPACE::close_memoryfile( mem_file );
            
            mem_file = ARCHIVE_WORKSPACE::open_memoryfile(identifier.c_str(), O_WRONLY);
            data = ARCHIVE_WORKSPACE::files(mem_file).data;
            ARCHIVE_WORKSPACE::close_memoryfile( mem_file );
            return;
        }
        else
            tar_skip_regfile(tar_ptr);
    }while(!th_read(tar_ptr));

    PHYSBAM_FATAL_ERROR("Can't read archive. File not found!");
      
}

void ARCHIVE::CloseArchive()
{
    //LOG::SCOPE scope("ARCHIVE::CloseArchive");

    if(!archive_ptr)
        PHYSBAM_FATAL_ERROR("Can't close archive. Archive not open!");

    TAR* tar_ptr = (TAR*)(archive_ptr);
    tar_close(tar_ptr);
}

void ARCHIVE::WriteArchive()
{
    //LOG::SCOPE scope("ARCHIVE::WriteArchive");

    if(!archive_ptr)
        PHYSBAM_FATAL_ERROR("Can't write archive. Archive not open!");

    std::ofstream tar_file_out(ARCHIVE_WORKSPACE::files(1).filename.c_str());
    tar_file_out.write( (char*)(ARCHIVE_WORKSPACE::files(1).data.Get_Array_Pointer()),
                        ARCHIVE_WORKSPACE::files(1).data.m);
    tar_file_out.close();
}

void ARCHIVE::CleanUpArchive()
{
    //LOG::SCOPE scope("ARCHIVE::CleanUpArchive");

    if(!archive_ptr)
        PHYSBAM_FATAL_ERROR("Can't cleanup archive. Archive not open!");

    archive_ptr = NULL;
    ARCHIVE_WORKSPACE::files.Resize(0);
}
