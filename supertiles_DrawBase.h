#ifndef __SUPERTILES_DRAW_BASE__
#define __SUPERTILES_DRAW_BASE__


template<typename DATA, typename DATA_OFF, typename DIMV>
struct supertiles_DrawBase
{
  using F = double;
  using F2 = V2<F>;
  
  auto getScale() const
  {
    return _scale;
  }

  auto getOffset() const
  {
    return _offset;
  }

  template<typename T>
  auto getTileDim(const T& nodeId) const
  {
    hassertm2(nodeId < _memberDims.size(), nodeId, _memberDims.size());
    return _memberDims[nodeId];
  }
  
  const DATA _data;
  const DATA_OFF _dataOffsets;
  const DIMV _memberDims;
  const F2 _offset;
  const F _scale;
  const uint32_t _nLOD;

protected:
  
  supertiles_DrawBase(const DATA data,
		      const DATA_OFF dataOffsets,
		      const DIMV memberDims,
		      const F2 offset, const F scale,
		      const uint32_t nLOD=1) :
    _data(data),
    _dataOffsets(dataOffsets),
    _memberDims(memberDims),
    _offset(offset),
    _scale(scale),
    _nLOD(nLOD)
  {
    assert(!_dataOffsets.empty());
    hassertm2(_dataOffsets.back()==data.size(), _dataOffsets.back(), data.size());        
  }
  
};

#endif //__SUPERTILES_DRAW_BASE__
