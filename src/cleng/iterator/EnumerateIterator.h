#include <iterator>

template<typename IteratorT>
class EnumerateIterator : std::iterator<std::forward_iterator_tag,
        typename std::iterator_traits<IteratorT>::value_type> {

private:
    size_t mCurIdx;
    IteratorT mItr;

public:
    using ValueT = typename std::iterator_traits<IteratorT>::value_type;
    using IdxValPair = std::pair<size_t, ValueT>;

    explicit EnumerateIterator(IteratorT &&iterator)
            : mCurIdx{0}, mItr{std::forward<IteratorT>(iterator)} {}

    EnumerateIterator(IteratorT &&iterator, size_t startingCount)
            : mCurIdx{startingCount}, mItr{std::forward<IteratorT>(iterator)} {}

    EnumerateIterator &operator++() {
        ++mItr;
        ++mCurIdx;
        return *this;
    }

    EnumerateIterator operator++(int) {
        auto temp{*this};
        operator++();
        return temp;
    }

    bool operator==(const EnumerateIterator &enumItr) const {
        return (mCurIdx == enumItr.mCurIdx) && (mItr == enumItr.mItr);
    }

    bool operator!=(const EnumerateIterator &enumItr) const {
        return !(*this == enumItr);
    }

    IdxValPair operator*() {
        return std::make_pair(mCurIdx, *mItr);
    }

    IdxValPair operator*() const {
        return std::make_pair(mCurIdx, *mItr);
    }

};

template<typename T>
struct EnumerateWrapper {
    T &range;
};

template<typename T>
EnumerateWrapper<T> Enumerate(T &&range) { return {range}; }

template<typename T>
auto begin(EnumerateWrapper<T> wrapper) {
    return EnumerateIterator<decltype(std::begin(wrapper.range))>(
            std::begin(wrapper.range));
}

template<typename T>
auto end(EnumerateWrapper<T> wrapper) {
    return EnumerateIterator<decltype(std::end(wrapper.range))>(
            std::end(wrapper.range),
            std::distance(std::begin(wrapper.range),
                          std::end(wrapper.range)));
}
