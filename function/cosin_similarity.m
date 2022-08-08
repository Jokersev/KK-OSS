function similarity = cosin_similarity(image1,image2)
% similarity = sum(dot(image1,image2))/(norm(image1)*norm(image2));
similarity = sum(dot(image1,image2))/(sqrt(sum(sum(image1.^2)))*sqrt(sum(sum(image2.^2))));
end

