#ifndef LOADFILE_H
#define LOADFILE_H

#include <vector>
#include <string>
#include "loadfiles_export.h"

class LOADFILES_EXPORT LoadFile
{
private:
    std::string sourcePath, targetPath, brand, model;
    struct FileInfo
    {
        std::string jsonPath;
        std::string filePath;
        std::string bone, side, wedge, implant_type, implant_name, json_name;
        bool exist;
        bool ExistPath(const std::string& pBrand, const std::string& pModel, const std::string& pDir) const;
        void saveFile(const std::string& pBrand, const std::string& pModel, const std::string& pDir, bool pOverwrite, std::vector<std::string>& errorList) const;
        std::string getNewJsonName() const;
    };
    std::vector<FileInfo> fileInfoList;

public:
    LoadFile(const std::string& pBrand, const std::string& pModel, const std::string& pSourcePath, const std::string& pTargetPath);
    
    /*
        This function reads and loads the files in the source folder into memory. 
        If a json file has a problem in its data then it is not loaded into memory 
        and the error is reported in errorList parameter.
        If any json file within the origin source folder exists in the target folder,
        it is reported in overwriteWarning parameter.
    */
    void readFiles(std::vector<std::string>& errorList, std::vector<std::string>& overwriteWarning);

    /*
        This function save the files loaded in memory in the source folder.
        If the save process have any problem is reported in errorList parameter.
        The pOverwrite parameter specifies whether to overwrite existing files or not,
        by default is true.
    */
    void writeFiles(std::vector<std::string>& errorList, bool pOverwrite = true) const;
};

#endif