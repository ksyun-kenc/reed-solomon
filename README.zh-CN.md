# reed-solomon

Reed-Solomon Erasure Coding in C++ used by Regame

## 目前版本说明

1. 支持 Cauchy、ParV1、Vandermonde 三种矩阵编码。
2. 实现 Encode、Reconstruct 函数。Verify 函数暂没使用，所以没实现，其实现和 Encode 类似，只是为了在有校验失败时尽快返回，所以数据处理顺序可能不同。
3. Encode 函数在计算矩阵时，先行后列和先列后行，性能都是不一样的，最多可能有 6 种组合，目前只实现普遍认为性能比较好的一种。
4. 提高可读性的做法：GF 的加减法都是 xor，但数学上到底是加还是减，意义不同，我的代码里用封装区分了，其它实现都是直接 xor。
5. 矩阵求逆有 Gaussian Elimination 和 Gauss-Jordan Elimination 两种算法，我用的是前者，有些实现并没有写明用那种算法。

## 后续优化方向

1. GF 运算目前是 C++ 实现的，可以调用 GF Complete，它有使用 SIMD 加速：<https://github.com/imp/gf-complete>
2. 矩阵运算是可以并行的，每行列运算都是独立的。
3. 矩阵求逆代价高，所以可以做一下关键矩阵的缓存层。
