# Integration of OpenMS library into your own project

This is the configuration that worked for me

## Building and Installing OpenMS

Following https://abibuilder.informatik.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/install_linux.html.
The installation that ended up working, was to install all the necessary tools via linux apt and not building them from contrib-lib.
Although it might be necessary if some packages do not work properly and I did a lot of reinstalling before it ended up working.

I also add the flag *-fPIC* to the CMakeList Compiler flags (unsure if necessary).
* set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CF_OPENMS_ADDCXX_FLAGS} -fPIC")


Then Building in an directory OpenMS-build (next to OpenMS dir) without linking contrib-build:
* cmake -DBOOST_USE_STATIC=OFF ../OpenMS
* make OpenMS
* make
* make test

These make commands build the bases for succesful installation. First only OpenMS build worked. Later make (all) did as well. Only then the test can be made.
99% Passed for me.

Then add
* export LD_LIBRARY_PATH="/home/ynowatzk/repos/OpenMS-build/lib:$LD_LIBRARY_PATH"

into .basrc in home directory.




## Including OpenMS

Make sure OpenMS is installed properly (see previous step). Then use the CMakeList.txt from the example external project within OpenMS directory.
Here the main issues came from finding the libraries Qt5 and OpenMP.
What made it work eventually was adding:
* find_package(Qt5 COMPONENTS Core Network REQUIRED)
* find_package(OpenMP)
to the CMakeList.txt.

Additionally, I explicitly linked the Path my OpenMS-build installation:
* find_package(OpenMS PATHS "/home/ynowatzk/repos/OpenMS-build")

And had to add *-fPIC* flag:
* set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMS_ADDCXX_FLAGS} -fPIC")

At last this cmake command ended up building successfully:
* cmake -D OpenMS_DIR=/home/ynowatzk/repos/OpenMS-build/ /home/ynowatzk/repos/OpenMS/share/OpenMS/examples/external_code
* make
* ./Main

Note: To get it running in the Clion IDE I had to add : -D OpenMS_DIR=/home/ynowatzk/repos/OpenMS-build/ to the Build,Execution,Deployment->Cmake->Cmake Options in the settings.
If all goes well, it Main says: All good and well!
