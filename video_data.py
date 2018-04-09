from fenics import Expression
import numpy as np
import imageio


def rgb_to_grayscale(rgbimage):
    return np.mean(rgbimage, axis=2)


class VideoData(Expression):
    def load_video(self, filename):
        self.t = 0
        self.frame = 0
        self.loop = True

        self.filename = filename
        self.video = imageio.get_reader(filename, "ffmpeg")
        meta = self.video.get_meta_data()
        self.fps = meta["fps"]
        self.dt = 1. / self.fps

        self.nframes = meta["nframes"]
        self.t_max = self.dt * self.nframes

        # IMPORTANT:  We change x and y coordinates, to have
        # right-handed coordinates
        self.shape = (meta["size"][1],meta["size"][0])

        self.current_image = rgb_to_grayscale(self.video.get_data(0))

        self.ref_point = np.array([0, 0])
        self.size = np.array([1, 1])

        self.last_frames = dict()

    def __str__(self):
        text = ""
        if self.video is not None:
            text = "\nVideo: %s, Length: %fs, Resolution: (%d, %d)" % (self.filename,
                                                                    self.t_max,
                                                                    *self.shape)
        return text

    def value_shape(self):
        return (1)

    def set_reference_frame(self, left_bottom_corner, size):
        self.ref_point = np.array(left_bottom_corner)
        self.size = np.array(size)

    def set_time(self, t):
        assert self.video is not None

        if t > self.t + self.dt or t < self.t:

            #if t < self.t:
            #    print("waring: read video frames in reverse!")

            i = int(np.floor(t/self.dt))

            if self.loop:
                i = i % (2*self.nframes)
                if i >= self.nframes:
                    i = 2*self.nframes - i

            if 0 <= i < self.nframes:
                self.frame = i

                if i in self.last_frames:
                    self.current_image = self.last_frames[i]
                else:
                    self.current_image = rgb_to_grayscale(self.video.get_data(self.frame))
                    self.last_frames[i] = self.current_image





    def point_to_index(self, x):
        assert self.video is not None

        # warning: only supports 2D!
        #
        # compute coordinates of the point with respect to
        # a coordinate system with the image spanning the box between (0,0) - (1,1)
        x = np.array(x)
        x = x - self.ref_point
        y = x / self.size
        # note that for a image the origin lies in the top left corner
        y[[1, 0]] = y[[0,1]]
        y[0] = 1. - y[0]

        # convert to indices
        index = y * np.array(self.shape)
        index = np.floor(index)
        return int(index[0]), int(index[1])

    def eval(self, value, x):
        assert self.current_image is not None

        ix, iy = self.point_to_index(x)

        if ix == self.shape[0]: ix -= 1
        if iy == self.shape[1]: iy -= 1

        if 0 <= ix < self.shape[0] and 0 <= iy < self.shape[1]:
            value[0] = self.current_image[ix,iy]
        else:
            print("Accessing index out of range (%d,%d)" % (ix,iy))
            value[0] = 0.
