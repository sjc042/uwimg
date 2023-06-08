from uwimg import *
# im = load_image("data/dog.jpg")
# f = make_gaussian_filter(2)
# lfreq = convolve_image(im, f, 1)
# hfreq = im - lfreq
# reconstruct = lfreq + hfreq
# print(f"og img [{im.c}x{im.h}x{im.w}]\n")
# print(f"low freq img [{lfreq.c}x{lfreq.h}x{lfreq.w}]\n")
# print(f"high freq img [{hfreq.c}x{hfreq.h}x{hfreq.w}]\n")
# print(f"reconstruct img [{reconstruct.c}x{reconstruct.h}x{reconstruct.w}]\n")
# save_image(lfreq, "low-frequency")
# save_image(hfreq, "high-frequency")
# save_image(reconstruct, "reconstruct")

im = load_image("data/dog.jpg")
# res = sobel_image(im)
# mag = res[0]
# feature_normalize(mag)
# save_image(mag, "magnitude")
im_sobelColor = colorize_sobel(im)
save_image(im_sobelColor, "im_sobelColor")