#include "LoadFiles.hpp"
#include <qfile.h>
#include <qdir.h>
#include <qstring.h>
#include <qjsondocument.h>
#include <qjsonobject.h>
#include <qjsonvalue.h>
#include <iostream>

inline bool checkFile(QDir& directory, const QString& filename)
{
    if (directory.exists(filename))
    {
        QString fullPath = directory.filePath(filename);
        directory.setPath(fullPath);
        return true;
    }
    else
    {
        return false;
    }
}

inline void checkAndCreateDir(QDir& directory, const QString& filename)
{
    if (directory.exists(filename) == false)
    {
        directory.mkdir(filename);
    }

    QString fullPath = directory.filePath(filename);
    directory.setPath(fullPath);
}

LoadFile::LoadFile(const std::string& pBrand, const std::string& pModel, const std::string& pSourcePath, const std::string& pTargetPath)
{
    sourcePath = pSourcePath;
    targetPath = pTargetPath;
    brand = pBrand;
    model = pModel;
}

void LoadFile::writeFiles(std::vector<std::string>& errorList, bool pOverwrite) const
{
    auto it1 = fileInfoList.begin();
    auto it2 = fileInfoList.end();

    for (; it1 != it2; ++it1)
    {
        it1->saveFile(brand, model, targetPath, pOverwrite, errorList);
    }
}

void LoadFile::readFiles(std::vector<std::string>& errorList, std::vector<std::string>& overwriteWarning)
{
    fileInfoList.clear();
    QDir directory(QString::fromStdString(sourcePath));
    QStringList files = directory.entryList(QStringList() << "*.json" << "*.JSON", QDir::Files);
    for each (QString filename in files)
    {
        QString jsonPath = directory.filePath(filename);
        QString var;
        QFile file;
        file.setFileName(jsonPath);
        file.open(QIODevice::ReadOnly | QIODevice::Text);
        var = file.readAll();
        file.close();
        QJsonDocument jsonDoc = QJsonDocument::fromJson(var.toUtf8());
        QJsonObject jsonObject = jsonDoc.object();

        QJsonValue jsonStl = jsonObject.value(QString("filename"));

        if (jsonStl.isNull())
        {
            std::string msg = "The " + filename.toStdString() + " file does not have 3D model file. Check the 'filename' field.";
            errorList.push_back(msg);
            continue;
        }

        if (directory.exists(jsonStl.toString()))
        {
            FileInfo myObject;
            QJsonValue jsonBone = jsonObject.value(QString("bone"));
            QJsonValue jsonType = jsonObject.value(QString("implant_type"));
            QJsonValue jsonSide;

            if (jsonBone.isNull())
            {
                std::string msg = "The " + filename.toStdString() + " file does not specify the type of bone. Check the 'bone' field.";
                errorList.push_back(msg);
                continue;
            }

            if (jsonType.isNull())
            {
                std::string msg = "The " + filename.toStdString() + " file does not specify the implant type. Check the 'implant_type' field.";
                errorList.push_back(msg);
                continue;
            }

            if (jsonBone.toString() != "patella")
            {
                jsonSide = jsonObject.value(QString("side"));

                if (jsonSide.isNull())
                {
                    std::string msg = "The " + filename.toStdString() + " file does not specify the bone side. Check the 'side' field.";
                    errorList.push_back(msg);
                    continue;
                }
            }

            QString implantType;

            if (jsonType.isString())
            {
                implantType = jsonType.toString();
            }
            else if (jsonType.isDouble())
            {
                int num = jsonType.toDouble();
                implantType = QString::number(num);
            }
            else
            {
                std::string msg = "The " + filename.toStdString() + " file does not have a right implant type format. Check the 'implant_type' field.";
                errorList.push_back(msg);
                continue;
            }

            QString stlPath = directory.filePath(jsonStl.toString());
            QString bone = jsonBone.toString();
            QString side;
            
            if (bone != "patella")
            {
                side = jsonSide.toString();
            }

            if (bone == "femur")
            {
                myObject.jsonPath = jsonPath.toStdString();
                myObject.filePath = stlPath.toStdString();
                myObject.bone = bone.toStdString();
                myObject.side = side.toStdString();
                myObject.implant_type = implantType.toStdString();
            }
            else if (bone == "patella")
            {
                myObject.jsonPath = jsonPath.toStdString();
                myObject.filePath = stlPath.toStdString();
                myObject.bone = bone.toStdString();
                myObject.implant_type = implantType.toStdString();
            }
            else if (bone == "tibia")
            {
                QJsonValue jsonWedge = jsonObject.value(QString("wedge"));
                if (jsonWedge.isNull())
                {
                    std::string msg = "The " + filename.toStdString() + " file does not specify the wedge size. Check the 'wedge' field.";
                    errorList.push_back(msg);
                    continue;
                }

                QString wedge;
                
                if (jsonWedge.isString())
                {
                    wedge = jsonWedge.toString();
                }
                else if (jsonWedge.isDouble())
                {
                    int num = jsonWedge.toDouble();
                    wedge = QString::number(num);
                }
                else
                {
                    std::string msg = "The " + filename.toStdString() + " file does not have a right wedge size format. Check the 'wedge' field.";
                    errorList.push_back(msg);
                    continue;
                }
                
                myObject.jsonPath = jsonPath.toStdString();
                myObject.filePath = stlPath.toStdString();
                myObject.bone = bone.toStdString();
                myObject.side = side.toStdString();
                myObject.wedge = wedge.toStdString();
                myObject.implant_type = implantType.toStdString();
            }
            else
            {
                std::string msg = "The " + filename.toStdString() + " file does not specify the type of bone. Check the 'bone' field.";
                errorList.push_back(msg);
                continue;
            }

            myObject.implant_name = jsonStl.toString().toStdString();
            myObject.json_name = filename.toStdString();
            myObject.exist = myObject.ExistPath(brand, model, targetPath);

            if (myObject.exist == true)
            {
                std::string msg = "The data in " + filename.toStdString() + " file already exists. Do you want to overwrite it?";
                overwriteWarning.push_back(msg);
            }

            fileInfoList.push_back(myObject);
        }
        else
        {
            QString wrongFile = jsonStl.toString();
            std::string msg = "The " + wrongFile.toStdString() + " file included in " + filename.toStdString() + " data file does not exist.";
            errorList.push_back(msg);
        }
    }
}

std::string LoadFile::FileInfo::getNewJsonName() const
{
    QString newJsonName;
    if (bone == "femur")
    {
        newJsonName = "femur_" + QString::fromStdString(side) + "_" + QString::fromStdString(implant_type) + "_data.json";
    }
    else if (bone == "patella")
    {
        newJsonName = "patella_" + QString::fromStdString(implant_type) + "_data.json";
    }
    else
    {
        newJsonName = "tibia_" + QString::fromStdString(implant_type) + "_" + QString::fromStdString(wedge) + "_data.json";
    }

    return newJsonName.toStdString();
}


void LoadFile::FileInfo::saveFile(const std::string& pBrand, const std::string& pModel, const std::string& pDir, bool pOverwrite, std::vector<std::string>& errorList) const
{
    if (pOverwrite == true || (pOverwrite == false && exist == false))
    {
        QDir directory(QString::fromStdString(pDir));
        checkAndCreateDir(directory, QString::fromStdString(pBrand));
        checkAndCreateDir(directory, QString::fromStdString(pModel));
        checkAndCreateDir(directory, QString::fromStdString(implant_type));
        checkAndCreateDir(directory, QString::fromStdString(bone));

        if (bone == "femur")
        {
            checkAndCreateDir(directory, QString::fromStdString(side));
        }
        else if (bone == "tibia")
        {
            checkAndCreateDir(directory, QString::fromStdString(side));
            checkAndCreateDir(directory, QString::fromStdString(wedge));
        }

        QStringList files = directory.entryList(QStringList() << "*.*", QDir::Files);
        for each (QString filename in files)
        {
            directory.remove(filename);
        }
        
        QString newJsonName = QString::fromStdString(getNewJsonName());
        QString newFileName = QString::fromStdString(implant_name);

        QString targetJsonPath = directory.filePath(newJsonName);
        QString targetFilePath = directory.filePath(newFileName);

        if (QFile::copy(QString::fromStdString(jsonPath), targetJsonPath) == false)
        {
            std::string msg = "Error copying the " + json_name + " file";
            errorList.push_back(msg);
        }

        if (QFile::copy(QString::fromStdString(filePath), targetFilePath) == false)
        {
            std::string msg = "Error copying the " + implant_name + " file";
            errorList.push_back(msg);
        }
    }
}

bool LoadFile::FileInfo::ExistPath(const std::string& pBrand, const std::string& pModel, const std::string& pDir) const
{
    QDir directory(QString::fromStdString(pDir));

    QString myBrand = QString::fromStdString(pBrand);
    QString myModel = QString::fromStdString(pModel);

    if (checkFile(directory, myBrand) == false)
    {
        return false;
    }

    if (checkFile(directory, myModel) == false)
    {
        return false;
    }

    if (checkFile(directory, QString::fromStdString(implant_type)) == false)
    {
        return false;
    }

    if (checkFile(directory, QString::fromStdString(bone)) == false)
    {
        return false;
    }

    if (bone == "femur")
    {
        if (checkFile(directory, QString::fromStdString(side)) == false)
        {
            return false;
        }

        QStringList files = directory.entryList(QStringList() << "*.json" << "*.JSON", QDir::Files);
        if (files.size() == 0)
        {
            return false;
        }
    }
    else
    {
        if (checkFile(directory, QString::fromStdString(wedge)) == false)
        {
            return false;
        }

        QStringList files = directory.entryList(QStringList() << "*.json" << "*.JSON", QDir::Files);
        if (files.size() == 0)
        {
            return false;
        }
    }

    return true;
}

