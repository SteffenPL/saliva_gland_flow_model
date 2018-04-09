from fenics import Mesh, Expression, FunctionSpace, plot, UnitSquareMesh, File

import matplotlib.pyplot as plt

from video_data import VideoData


def test_video_loading():
    m = UnitSquareMesh(30,30)
    V = FunctionSpace(m, "CG", 4)

    video = VideoData(element=V.ufl_element())
    video.load_video("video_data/ach12.mp4")

    print( video )

    plot(video,mesh=m)
    plt.show()

    video.set_time(10000.)
    plot(video, mesh=m)
    plt.show()

    assert True