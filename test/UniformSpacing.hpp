
    template<typename ImageType>
    bool UniformSpacing(const std::vector<std::string>& pFiles)
    {
        auto it1 = pFiles.begin();
        auto it2 = pFiles.end();
        float previous, current, distPrev = -1, distCurrent = -1, thickness = 0;
        bool sameThickness = true;
        for (; it1 != it2; ++it1)
        {
            using ReaderImage = itk::ImageFileReader<ImageType>;
            typename ReaderImage::Pointer reader = ReaderImage::New();
            reader->SetFileName(*it1);
            reader->Update();
            ImageType::Pointer image = reader->GetOutput();
            ImageType::PointType index = image->GetOrigin();

            if (it1 == pFiles.begin())
            {
                previous = index[2];
                current = index[2];
            }
            else
            {
                current = index[2];
                distCurrent = abs(current - previous);
                previous = current;

                if (distPrev < 0)
                {
                    distPrev = distCurrent;
                }

                thickness = abs(distCurrent - distPrev);
                distPrev = distCurrent;

                if (thickness > 0.05)
                {
                    sameThickness = false;
                    break;
                }
            }
        }

        return sameThickness;
    }
