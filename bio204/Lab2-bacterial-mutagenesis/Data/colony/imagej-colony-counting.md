# 使用ImageJ数菌落

1. 打开图片, 注意数菌落的图片要**裁掉培养皿边缘**, 不然会误判
2. `Image` > `Type` > `8-bit`
3. `Image` > `Adjust` > `Brightness/Contrast...` > `Auto` > `Apply` > `OK`
4. `Edit` > `Invert`
5.  `Image` > `Adjust` > `Threshold...` > 设定范围, 尽量使菌落为红色, 其他均不为红色, 尤其是培养皿边缘, 且菌落见尽量不连接 > `Apply` 
6. `Process` > `Binary` > `Watershed`
7. `Analyze`  > `Analyze Particles...` > 设定尺寸和圆度范围, Show: `Outlines` > `OK`
8. 在`ROI Manager`中调整, `Flatten [F]`渲染所有图层, 右键后`Control Panel...`保存