def orb_align_images(ref_img, test_img, n_features=4000, draw_matches=True):
    ref = cv2.cvtColor(ref_img, cv2.COLOR_BGR2GRAY)
    test = cv2.cvtColor(test_img, cv2.COLOR_BGR2GRAY)

    orb = cv2.ORB_create(n_features)
    kp1, des1 = orb.detectAndCompute(ref, None)
    kp2, des2 = orb.detectAndCompute(test, None)

    matcher = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=True)
    matches = matcher.match(des1, des2)
    matches = sorted(matches, key=lambda x: x.distance)

    src_pts = np.float32([kp1[m.queryIdx].pt for m in matches]).reshape(-1,1,2)
    dst_pts = np.float32([kp2[m.trainIdx].pt for m in matches]).reshape(-1,1,2)

    H, mask = cv2.findHomography(dst_pts, src_pts, cv2.RANSAC, 5.0)

    aligned_color = cv2.warpPerspective(test_img, H, (ref_img.shape[1], ref_img.shape[0]))

    if draw_matches:
        vis = cv2.drawMatches(ref, kp1, test, kp2, matches[:50], None, flags=2)
        plt.figure(figsize=(16,7))
        plt.subplot(1,2,1); plt.imshow(vis, cmap='gray')
        plt.subplot(1,2,2); plt.imshow(cv2.cvtColor(aligned_color, cv2.COLOR_BGR2RGB))
        plt.show()

    return aligned_color, H