from uwimg import *

# a = load_image("data/dog_a.jpg")
# b = load_image("data/dog_b.jpg")
# flow = optical_flow_images(b, a, 15, 8)
# draw_flow(a, flow, 8)
# save_image(a, "lines")

# optical_flow_webcam(15,4,8)

import cv2
cap = cv2.VideoCapture(0)
# cap2 = cv2.VideoCapture(0)
while True:

    ret, frame = cap.read()
    cv2.imshow('frame',frame)
    if cv2.waitKey(1) & 0xFF == ord('q'):
        break

cap.release()
cv2.destroyAllWindows()
