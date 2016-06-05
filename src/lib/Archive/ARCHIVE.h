#include <PhysBAM_Tools/Arrays/ARRAY.h>

namespace PhysBAM{

    typedef ARRAY<char> DATA_BUFFER;

    struct ARCHIVE
    {
        static void CreateArchive(std::string archive_name);
        static void OpenArchive(std::string archive_name);
        
        static void WriteDataBuffer(std::string identifier, const DATA_BUFFER& data);
        static void ReadDataBuffer(std::string identifier, DATA_BUFFER& data);
        
        static void CloseArchive();
        static void WriteArchive();
        static void CleanUpArchive();

    private:
        
        static void* archive_ptr;
    };

};
