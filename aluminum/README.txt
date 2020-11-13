sweep-3: 加入了surface reaction

sweep-2: 加入inception后，不能运行，发生CVODE ERROR。
但是，去掉XML文件中的sintering后，便能顺利运行。

已做过将kinetics中的A值减小的试验，结果良好。
但，A的值本身是对的，不需要减少。
A大，说明rate就相应地大，rxn就快。
将时间间隔无限减小有效果吗。


chem-diffusion: 利用Arrhenius equation代表diffusion。
